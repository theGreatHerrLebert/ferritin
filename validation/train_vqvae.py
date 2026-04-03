#!/usr/bin/env python3
"""Train a VQ-VAE structural alphabet encoder from scratch.

No FoldSeek dependency — the 20 structural states are learned purely from
backbone geometry via vector quantization, matching the approach described
in van Kempen et al. Nature Biotechnology (2023).

Architecture:
  Encoder: 10 → 10 (ReLU) → 10 (ReLU) → 2 (linear)
  VQ: 20 codebook vectors in 2D
  Decoder: 2 → 10 (ReLU) → 10 (ReLU) → 10 (linear)

Loss: reconstruction_loss + β * commitment_loss

The 20 states emerge from the VQ bottleneck — no labels needed.

Usage:
    python validation/train_vqvae.py --pdb-dir validation/pdbs_10k/ --n-structures 10000
    python validation/train_vqvae.py --pdb-dir validation/pdbs_10k/ --n-structures 10000 --epochs 1000
"""

import argparse
import json
import os
import random
import sys
import time

import numpy as np

# ============================================================================
# Feature extraction (pure Python, matches our Rust implementation)
# ============================================================================

DISTANCE_CA_CB = 1.5336


def vnorm(a):
    n = np.linalg.norm(a)
    return a / n if n > 1e-12 else np.zeros(3)


def approx_cb(ca, n, c):
    v1 = vnorm(c - ca)
    v2 = vnorm(n - ca)
    b1 = v2 + v1 / 3.0
    b2 = np.cross(v1, b1)
    u1 = vnorm(b1)
    u2 = vnorm(b2)
    term1 = -v1 / 3.0
    term2a = -u1 / 2.0 - u2 * np.sqrt(3) / 2.0
    term2 = term2a * np.sqrt(8) / 3.0
    v4 = term1 + term2
    return ca + v4 * DISTANCE_CA_CB


def rodrigues_rotate(v, k, angle):
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    return v * cos_a + np.cross(k, v) * sin_a + k * np.dot(k, v) * (1 - cos_a)


def virtual_center(ca, cb, n):
    alpha = np.radians(270.0)
    beta = np.radians(0.0)
    d = 2.0
    v = cb - ca
    a = cb - ca
    b = n - ca
    k = vnorm(np.cross(a, b))
    v = rodrigues_rotate(v, k, alpha)
    k2 = vnorm(n - ca)
    v = rodrigues_rotate(v, k2, beta)
    return ca + v * d


def compute_features(ca, i, j):
    u1 = vnorm(ca[i] - ca[i - 1])
    u2 = vnorm(ca[i + 1] - ca[i])
    u3 = vnorm(ca[j] - ca[j - 1])
    u4 = vnorm(ca[j + 1] - ca[j])
    u5 = vnorm(ca[j] - ca[i])
    sep = float(j - i)
    return [
        np.dot(u1, u2),
        np.dot(u3, u4),
        np.dot(u1, u5),
        np.dot(u3, u5),
        np.dot(u1, u4),
        np.dot(u2, u3),
        np.dot(u1, u3),
        np.linalg.norm(ca[i] - ca[j]),
        np.sign(sep) * min(abs(sep), 4.0),
        np.sign(sep) * np.log(abs(sep) + 1.0),
    ]


def extract_backbone(pdb_path):
    """Extract CA, N, C, CB coordinates from a PDB file."""
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name not in ("CA", "N", "C", "CB"):
                continue
            chain = line[21]
            resnum = line[22:27].strip()
            key = (chain, resnum)
            if key not in residues:
                residues[key] = {}
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            residues[key][atom_name] = np.array([x, y, z])

    ca, n_atoms, c_atoms, cb_atoms = [], [], [], []
    for key in sorted(residues.keys()):
        atoms = residues[key]
        if "CA" not in atoms or "N" not in atoms or "C" not in atoms:
            continue
        ca.append(atoms["CA"])
        n_atoms.append(atoms["N"])
        c_atoms.append(atoms["C"])
        if "CB" in atoms:
            cb_atoms.append(atoms["CB"])
        else:
            cb_atoms.append(approx_cb(atoms["CA"], atoms["N"], atoms["C"]))

    if len(ca) < 3:
        return None
    return np.array(ca), np.array(n_atoms), np.array(c_atoms), np.array(cb_atoms)


def compute_all_features(ca, n_atoms, c_atoms, cb_atoms):
    """Compute features for all residues."""
    length = len(ca)
    virt = np.array([virtual_center(ca[i], cb_atoms[i], n_atoms[i]) for i in range(length)])

    partners = [-1] * length
    for i in range(1, length - 1):
        dists = np.linalg.norm(virt[1:length-1] - virt[i], axis=1)
        dists[i - 1] = float("inf")
        partners[i] = int(np.argmin(dists)) + 1

    features = []
    for i in range(1, length - 1):
        j = partners[i]
        if j < 1 or j >= length - 1:
            continue
        features.append(compute_features(ca, i, j))

    return np.array(features) if features else np.zeros((0, 10))


# ============================================================================
# VQ-VAE model (PyTorch)
# ============================================================================

def train_vqvae(features, n_states=20, embed_dim=2, epochs=500, beta=0.25):
    """Train a VQ-VAE on geometric features.

    The encoder maps 10D features → 2D embedding.
    VQ quantizes to nearest of 20 codebook vectors.
    The decoder reconstructs 10D features from the quantized embedding.

    Loss = ||x - x_hat||^2 + β * ||emb - sg(codebook)||^2
           + ||sg(emb) - codebook||^2

    No labels needed — the 20 states emerge from compression.
    """
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch.utils.data import DataLoader, TensorDataset

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"    Device: {device}")

    # Standardize
    mean = features.mean(axis=0)
    std = features.std(axis=0)
    std[std < 1e-8] = 1.0
    X = (features - mean) / std

    # Train/val split
    n = len(X)
    perm = np.random.RandomState(42).permutation(n)
    n_val = n // 10
    val_idx, train_idx = perm[:n_val], perm[n_val:]

    X_train = torch.tensor(X[train_idx], dtype=torch.float32, device=device)
    X_val = torch.tensor(X[val_idx], dtype=torch.float32, device=device)

    train_ds = TensorDataset(X_train)
    train_dl = DataLoader(train_ds, batch_size=1024, shuffle=True)

    # Encoder: 10 → 10 (ReLU) → 10 (ReLU) → embed_dim
    encoder = nn.Sequential(
        nn.Linear(10, 10), nn.ReLU(),
        nn.Linear(10, 10), nn.ReLU(),
        nn.Linear(10, embed_dim),
    ).to(device)

    # Decoder: embed_dim → 10 (ReLU) → 10 (ReLU) → 10
    decoder = nn.Sequential(
        nn.Linear(embed_dim, 10), nn.ReLU(),
        nn.Linear(10, 10), nn.ReLU(),
        nn.Linear(10, 10),
    ).to(device)

    # VQ codebook: 20 vectors in embed_dim space
    codebook = nn.Parameter(torch.randn(n_states, embed_dim, device=device) * 0.1)

    params = list(encoder.parameters()) + list(decoder.parameters()) + [codebook]
    optimizer = torch.optim.Adam(params, lr=1e-3)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    best_val_loss = float("inf")
    best_state = None

    for epoch in range(epochs):
        encoder.train(); decoder.train()
        epoch_recon = 0.0
        epoch_n = 0

        for (Xb,) in train_dl:
            # Encode
            emb = encoder(Xb)  # [B, 2]

            # VQ: find nearest codebook vector
            dists = torch.cdist(emb, codebook)  # [B, 20]
            indices = dists.argmin(dim=1)        # [B]
            quantized = codebook[indices]        # [B, 2]

            # Straight-through estimator: gradient flows through quantized
            emb_st = emb + (quantized - emb).detach()

            # Decode
            x_hat = decoder(emb_st)  # [B, 10]

            # Losses
            recon_loss = F.mse_loss(x_hat, Xb)
            commitment_loss = F.mse_loss(emb, quantized.detach())
            codebook_loss = F.mse_loss(quantized, emb.detach())

            loss = recon_loss + beta * commitment_loss + codebook_loss

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            epoch_recon += recon_loss.item() * len(Xb)
            epoch_n += len(Xb)

        scheduler.step()

        if (epoch + 1) % 100 == 0 or epoch == 0:
            encoder.eval(); decoder.eval()
            with torch.no_grad():
                val_emb = encoder(X_val)
                val_dists = torch.cdist(val_emb, codebook)
                val_idx_q = val_dists.argmin(dim=1)
                val_q = codebook[val_idx_q]
                val_hat = decoder(val_q)
                val_recon = F.mse_loss(val_hat, X_val).item()

                # Codebook usage
                usage = len(set(val_idx_q.cpu().numpy().tolist()))

                # Perplexity (how uniformly codebook is used)
                counts = torch.bincount(val_idx_q, minlength=n_states).float()
                probs = counts / counts.sum()
                perplexity = torch.exp(-torch.sum(probs * torch.log(probs + 1e-10))).item()

            train_recon = epoch_recon / epoch_n
            print(f"    Epoch {epoch+1:4d}: train_recon={train_recon:.4f}  "
                  f"val_recon={val_recon:.4f}  usage={usage}/20  "
                  f"perplexity={perplexity:.1f}  lr={scheduler.get_last_lr()[0]:.1e}")

            if val_recon < best_val_loss:
                best_val_loss = val_recon
                best_state = {
                    "encoder": {k: v.cpu().clone() for k, v in encoder.state_dict().items()},
                    "decoder": {k: v.cpu().clone() for k, v in decoder.state_dict().items()},
                    "codebook": codebook.detach().cpu().clone(),
                }

        elif epoch % 20 == 0:
            encoder.eval(); decoder.eval()
            with torch.no_grad():
                val_emb = encoder(X_val)
                val_dists = torch.cdist(val_emb, codebook)
                val_idx_q = val_dists.argmin(dim=1)
                val_q = codebook[val_idx_q]
                val_hat = decoder(val_q)
                val_recon = F.mse_loss(val_hat, X_val).item()
            if val_recon < best_val_loss:
                best_val_loss = val_recon
                best_state = {
                    "encoder": {k: v.cpu().clone() for k, v in encoder.state_dict().items()},
                    "decoder": {k: v.cpu().clone() for k, v in decoder.state_dict().items()},
                    "codebook": codebook.detach().cpu().clone(),
                }

    # Restore best
    encoder.load_state_dict(best_state["encoder"])
    encoder.eval(); encoder.cpu()
    best_codebook = best_state["codebook"].numpy()

    # Extract encoder weights
    state = best_state["encoder"]
    W1 = state["0.weight"].numpy()
    b1 = state["0.bias"].numpy()
    W2 = state["2.weight"].numpy()
    b2 = state["2.bias"].numpy()
    W3 = state["4.weight"].numpy()
    b3 = state["4.bias"].numpy()

    # Fold standardization into first layer
    W1_eff = W1 / std[np.newaxis, :]
    b1_eff = b1 - W1_eff @ mean

    # Final stats
    X_all = torch.tensor(X, dtype=torch.float32)
    with torch.no_grad():
        all_emb = encoder(X_all).numpy()
    dists = np.sum((all_emb[:, np.newaxis, :] - best_codebook[np.newaxis, :, :]) ** 2, axis=2)
    all_states = dists.argmin(axis=1)
    usage = len(set(all_states.tolist()))
    counts = np.bincount(all_states, minlength=n_states)
    probs = counts / counts.sum()
    perplexity = np.exp(-np.sum(probs * np.log(probs + 1e-10)))

    print(f"\n    Final: val_recon={best_val_loss:.4f}  "
          f"usage={usage}/20  perplexity={perplexity:.1f}")
    print(f"    State distribution: {dict(enumerate(counts.tolist()))}")

    return {
        "W1": W1_eff, "b1": b1_eff,
        "W2": W2, "b2": b2,
        "W3": W3, "b3": b3,
        "centroids": best_codebook,
        "val_recon_loss": float(best_val_loss),
        "codebook_usage": usage,
        "perplexity": float(perplexity),
        "mean": mean, "std": std,
        "embed_dim": embed_dim,
    }


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Train VQ-VAE structural alphabet")
    parser.add_argument("--pdb-dir", default="/scratch/TMAlign/test-pdbs/")
    parser.add_argument("--n-structures", type=int, default=100)
    parser.add_argument("--output", default="validation/alphabet_vqvae.json")
    parser.add_argument("--embed-dim", type=int, default=2)
    parser.add_argument("--epochs", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    all_pdbs = sorted([f for f in os.listdir(args.pdb_dir) if f.endswith(".pdb")])
    if args.n_structures < len(all_pdbs):
        random.seed(args.seed)
        pdbs = sorted(random.sample(all_pdbs, args.n_structures))
    else:
        pdbs = all_pdbs

    print(f"Training VQ-VAE on {len(pdbs)} structures")
    print(f"PDB dir: {args.pdb_dir}")
    print(f"Embed dim: {args.embed_dim}, Epochs: {args.epochs}")

    # Extract features from all structures
    print("\nExtracting geometric features...")
    all_features = []
    n_processed = 0
    n_skipped = 0
    t0 = time.time()

    for i, pdb_file in enumerate(pdbs):
        pdb_path = os.path.join(args.pdb_dir, pdb_file)
        try:
            backbone = extract_backbone(pdb_path)
        except Exception:
            n_skipped += 1
            continue

        if backbone is None:
            n_skipped += 1
            continue

        ca, n_atoms, c_atoms, cb_atoms = backbone
        features = compute_all_features(ca, n_atoms, c_atoms, cb_atoms)
        if len(features) > 0:
            all_features.append(features)
            n_processed += 1

        if (i + 1) % 500 == 0:
            elapsed = time.time() - t0
            n_res = sum(len(f) for f in all_features)
            print(f"  [{i+1}/{len(pdbs)}] {n_processed} structures, "
                  f"{n_res} residues, {elapsed:.0f}s")

    all_features = np.vstack(all_features)
    elapsed = time.time() - t0
    print(f"\n  {n_processed} structures, {len(all_features)} residues "
          f"({n_skipped} skipped) in {elapsed:.0f}s")
    print(f"  Feature shape: {all_features.shape}")

    if len(all_features) < 1000:
        print("ERROR: Not enough training data")
        sys.exit(1)

    # Train VQ-VAE
    print(f"\nTraining VQ-VAE (embed_dim={args.embed_dim}, epochs={args.epochs})...")
    trained = train_vqvae(all_features, embed_dim=args.embed_dim, epochs=args.epochs)

    # Save weights
    result = {
        "W1": trained["W1"].tolist(), "b1": trained["b1"].tolist(),
        "W2": trained["W2"].tolist(), "b2": trained["b2"].tolist(),
        "W3": trained["W3"].tolist(), "b3": trained["b3"].tolist(),
        "centroids": trained["centroids"].tolist(),
        "embed_dim": trained["embed_dim"],
        "val_recon_loss": trained["val_recon_loss"],
        "codebook_usage": trained["codebook_usage"],
        "perplexity": trained["perplexity"],
        "n_training_residues": int(len(all_features)),
        "n_structures": n_processed,
    }

    with open(args.output, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\n  Saved to {args.output}")

    # Generate Rust code
    W1, b1 = trained["W1"], trained["b1"]
    W2, b2 = trained["W2"], trained["b2"]
    W3, b3 = trained["W3"], trained["b3"]
    centroids = trained["centroids"]
    ed = trained["embed_dim"]

    rust_file = args.output.replace(".json", "_rust.txt")
    with open(rust_file, "w") as f:
        f.write("// Auto-generated VQ-VAE encoder weights for structural alphabet\n")
        f.write(f"// Trained on {n_processed} structures, "
                f"{len(all_features)} residues\n")
        f.write(f"// Architecture: 10→10(ReLU)→10(ReLU)→{ed}(linear)→VQ({20})\n")
        f.write(f"// Val recon loss: {trained['val_recon_loss']:.4f}, "
                f"perplexity: {trained['perplexity']:.1f}, "
                f"usage: {trained['codebook_usage']}/20\n")
        f.write(f"// Standardization folded into W1/b1\n\n")

        def write_matrix(name, mat):
            rows, cols = mat.shape
            f.write(f"static {name}: [[f64; {cols}]; {rows}] = [\n")
            for row in mat:
                vals = ", ".join(f"{v: .10f}" for v in row)
                f.write(f"    [{vals}],\n")
            f.write("];\n\n")

        def write_vector(name, vec):
            f.write(f"static {name}: [f64; {len(vec)}] = [\n    ")
            f.write(", ".join(f"{v: .10f}" for v in vec))
            f.write("\n];\n\n")

        write_matrix("ENCODER_W1", W1)
        write_vector("ENCODER_B1", b1)
        write_matrix("ENCODER_W2", W2)
        write_vector("ENCODER_B2", b2)
        write_matrix("ENCODER_W3", W3)
        write_vector("ENCODER_B3", b3)

        f.write(f"static CENTROIDS: [[f64; {ed}]; NUM_STATES] = [\n")
        for row in centroids:
            vals = ", ".join(f"{v: .10f}" for v in row)
            f.write(f"    [{vals}],\n")
        f.write("];\n")

    print(f"  Rust code saved to {rust_file}")


if __name__ == "__main__":
    main()
