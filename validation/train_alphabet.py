#!/usr/bin/env python3
"""Train structural alphabet encoder weights using FoldSeek as oracle.

Generates training data by:
1. Computing 10 geometric features from CA backbone coordinates
2. Running FoldSeek to get ground-truth 3Di state assignments
3. Training an MLP (PyTorch) to reproduce FoldSeek's 3Di output
4. Exporting weights for the Rust inference engine

The resulting weights are GPL-free since they are trained independently —
only the *algorithm* (feature extraction from backbone geometry) comes
from the published paper, not from FoldSeek's code.

Usage:
    python validation/train_alphabet.py --pdb-dir test-pdbs/ --n-structures 100
    python validation/train_alphabet.py --pdb-dir validation/pdbs/ --n-structures 1000
"""

import argparse
import json
import os
import random
import subprocess
import sys
import tempfile

import numpy as np

FOLDSEEK_BIN = "/scratch/TMAlign/foldseek/bin/foldseek"

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
    """Compute features + partners for all residues."""
    length = len(ca)
    virt = np.array([virtual_center(ca[i], cb_atoms[i], n_atoms[i]) for i in range(length)])

    partners = [-1] * length
    for i in range(1, length - 1):
        dists = np.linalg.norm(virt[1:length-1] - virt[i], axis=1)
        dists[i - 1] = float("inf")  # exclude self
        partners[i] = int(np.argmin(dists)) + 1

    features = []
    valid_indices = []
    for i in range(1, length - 1):
        j = partners[i]
        if j < 1 or j >= length - 1:
            continue
        features.append(compute_features(ca, i, j))
        valid_indices.append(i)

    return np.array(features) if features else np.zeros((0, 10)), valid_indices, partners


# ============================================================================
# FoldSeek oracle
# ============================================================================

def run_foldseek_3di(pdb_dir, pdb_files, tmpdir):
    """Run FoldSeek to get 3Di sequences for each structure."""
    dbdir = os.path.join(tmpdir, "db")
    result = subprocess.run(
        [FOLDSEEK_BIN, "createdb", pdb_dir, dbdir, "--threads", "4"],
        capture_output=True, text=True, timeout=120,
    )
    if result.returncode != 0:
        print(f"createdb failed: {result.stderr[:200]}")
        return {}

    subprocess.run(
        [FOLDSEEK_BIN, "lndb", dbdir + "_h", dbdir + "_ss_h"],
        capture_output=True, text=True, timeout=30,
    )

    outfile = os.path.join(tmpdir, "3di.fasta")
    result = subprocess.run(
        [FOLDSEEK_BIN, "convert2fasta", dbdir + "_ss", outfile],
        capture_output=True, text=True, timeout=60,
    )
    if result.returncode != 0:
        print(f"convert2fasta failed: {result.stderr[:200]}")
        return {}

    seqs = {}
    current_name = None
    current_seq = []
    with open(outfile) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name and current_seq:
                    seqs[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0].split(".")[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_name and current_seq:
        seqs[current_name] = "".join(current_seq)

    return seqs


# ============================================================================
# PyTorch training
# ============================================================================

def train_encoder(features, labels, n_states=20, embed_dim=2):
    """Train encoder matching FoldSeek arch: 10→10(ReLU)→10(ReLU)→2 + VQ.

    FoldSeek's actual architecture (decoded from kerasify binary):
      Dense(10→10, ReLU) → Dense(10→10, ReLU) → Dense(10→2, linear) → VQ(20)

    We train with cross-entropy + VQ commitment loss to learn an embedding
    space where the 20 structural states form separable clusters.
    """
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch.utils.data import DataLoader, TensorDataset

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"    Device: {device}")

    # Standardize features
    mean = features.mean(axis=0)
    std = features.std(axis=0)
    std[std < 1e-8] = 1.0
    X = (features - mean) / std

    # Train/val split
    n = len(X)
    perm = np.random.RandomState(42).permutation(n)
    n_val = n // 5
    val_idx, train_idx = perm[:n_val], perm[n_val:]

    X_train = torch.tensor(X[train_idx], dtype=torch.float32, device=device)
    y_train = torch.tensor(labels[train_idx], dtype=torch.long, device=device)
    X_val = torch.tensor(X[val_idx], dtype=torch.float32, device=device)
    y_val = torch.tensor(labels[val_idx], dtype=torch.long, device=device)

    train_ds = TensorDataset(X_train, y_train)
    train_dl = DataLoader(train_ds, batch_size=512, shuffle=True)

    # Direct classifier: 10 → 64 (ReLU) → 32 (ReLU) → 20
    # No 2D bottleneck — just classify directly. For Rust inference this is
    # 3 matrix multiplies + ReLU + argmax. Tiny and fast.
    model = nn.Sequential(
        nn.Linear(10, 64), nn.ReLU(), nn.Dropout(0.1),
        nn.Linear(64, 32), nn.ReLU(), nn.Dropout(0.1),
        nn.Linear(32, n_states),
    ).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=500)

    best_val_acc = 0.0
    best_state = None

    for epoch in range(500):
        model.train()
        for Xb, yb in train_dl:
            loss = F.cross_entropy(model(Xb), yb)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        scheduler.step()

        if (epoch + 1) % 100 == 0 or epoch == 0:
            model.eval()
            with torch.no_grad():
                val_acc = (model(X_val).argmax(1) == y_val).float().mean().item()
                train_acc = (model(X_train).argmax(1) == y_train).float().mean().item()
            print(f"    Epoch {epoch+1:3d}: train={train_acc:.3f}  val={val_acc:.3f}  "
                  f"lr={scheduler.get_last_lr()[0]:.1e}")
            if val_acc > best_val_acc:
                best_val_acc = val_acc
                best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
        elif epoch % 20 == 0:
            model.eval()
            with torch.no_grad():
                val_acc = (model(X_val).argmax(1) == y_val).float().mean().item()
            if val_acc > best_val_acc:
                best_val_acc = val_acc
                best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}

    # Restore best
    model.load_state_dict(best_state)
    model.eval(); model.cpu()

    # Extract weights (skip Dropout layers: indices 0,3,6)
    W1 = best_state["0.weight"].numpy()  # [64, 10]
    b1 = best_state["0.bias"].numpy()
    W2 = best_state["3.weight"].numpy()  # [32, 64]
    b2 = best_state["3.bias"].numpy()
    W3 = best_state["6.weight"].numpy()  # [20, 32]
    b3 = best_state["6.bias"].numpy()

    # Fold standardization into first layer
    W1_eff = W1 / std[np.newaxis, :]
    b1_eff = b1 - W1_eff @ mean

    # Final accuracy
    X_all = torch.tensor(X, dtype=torch.float32)
    with torch.no_grad():
        all_pred = model(X_all).argmax(1).numpy()
    accuracy = (all_pred == labels).mean()

    return {
        "W1": W1_eff, "b1": b1_eff,
        "W2": W2, "b2": b2,
        "W3": W3, "b3": b3,
        "centroids": None,  # no VQ — direct classification
        "accuracy": float(accuracy),
        "val_accuracy": float(best_val_acc),
        "mean": mean, "std": std,
        "embed_dim": n_states,  # output dim = 20
    }


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Train structural alphabet encoder")
    parser.add_argument("--pdb-dir", default="/scratch/TMAlign/test-pdbs/")
    parser.add_argument("--n-structures", type=int, default=100)
    parser.add_argument("--output", default="validation/alphabet_weights.json")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if not os.path.exists(FOLDSEEK_BIN):
        print(f"ERROR: FoldSeek not found at {FOLDSEEK_BIN}")
        sys.exit(1)

    all_pdbs = sorted([f for f in os.listdir(args.pdb_dir) if f.endswith(".pdb")])
    if args.n_structures < len(all_pdbs):
        random.seed(args.seed)
        pdbs = sorted(random.sample(all_pdbs, args.n_structures))
    else:
        pdbs = all_pdbs

    print(f"Training alphabet encoder on {len(pdbs)} structures")
    print(f"PDB dir: {args.pdb_dir}")

    # Step 1: Get FoldSeek 3Di assignments
    print("\nStep 1: Running FoldSeek to get ground-truth 3Di states...")
    with tempfile.TemporaryDirectory() as tmpdir:
        foldseek_seqs = run_foldseek_3di(args.pdb_dir, pdbs, tmpdir)
    print(f"  Got 3Di sequences for {len(foldseek_seqs)} structures")

    fs_alphabet = "ACDEFGHIKLMNPQRSTVWY"
    char_to_idx = {c: i for i, c in enumerate(fs_alphabet)}

    # Step 2: Extract features and match with FoldSeek labels
    print("\nStep 2: Extracting geometric features...")
    all_features = []
    all_labels = []
    n_matched = 0
    n_skipped = 0

    for pdb_file in pdbs:
        basename = os.path.splitext(pdb_file)[0]
        if basename not in foldseek_seqs:
            n_skipped += 1
            continue

        pdb_path = os.path.join(args.pdb_dir, pdb_file)
        backbone = extract_backbone(pdb_path)
        if backbone is None:
            n_skipped += 1
            continue

        ca, n_atoms, c_atoms, cb_atoms = backbone
        features, valid_indices, partners = compute_all_features(ca, n_atoms, c_atoms, cb_atoms)

        fs_seq = foldseek_seqs[basename]
        if len(fs_seq) != len(ca):
            n_skipped += 1
            continue

        for feat_idx, res_idx in enumerate(valid_indices):
            if res_idx >= len(fs_seq):
                continue
            fs_char = fs_seq[res_idx]
            if fs_char not in char_to_idx:
                continue
            all_features.append(features[feat_idx])
            all_labels.append(char_to_idx[fs_char])
            n_matched += 1

    all_features = np.array(all_features)
    all_labels = np.array(all_labels)

    print(f"  {n_matched} residues matched, {n_skipped} structures skipped")
    print(f"  Feature shape: {all_features.shape}")

    from collections import Counter
    dist = Counter(all_labels)
    print(f"  States represented: {len(dist)}/20")

    if len(all_features) < 100:
        print("ERROR: Not enough training data")
        sys.exit(1)

    # Step 3: Train
    print("\nStep 3: Training MLP encoder (PyTorch)...")
    trained = train_encoder(all_features, all_labels)

    accuracy = trained["accuracy"]
    val_accuracy = trained["val_accuracy"]
    print(f"\n  Final: train={accuracy:.1%}, val={val_accuracy:.1%}")

    # Step 4: Save weights as JSON
    result = {
        "W1": trained["W1"].tolist(), "b1": trained["b1"].tolist(),
        "W2": trained["W2"].tolist(), "b2": trained["b2"].tolist(),
        "W3": trained["W3"].tolist(), "b3": trained["b3"].tolist(),
        "centroids": trained["centroids"].tolist(),
        "embed_dim": trained["embed_dim"],
        "accuracy": round(accuracy, 4),
        "val_accuracy": round(val_accuracy, 4),
        "n_training_residues": int(len(all_features)),
        "n_structures": int(len(pdbs) - n_skipped),
        "alphabet": fs_alphabet,
    }

    with open(args.output, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  Saved weights to {args.output}")

    # Step 5: Generate Rust constants
    rust_file = args.output.replace(".json", "_rust.txt")
    W1, b1 = trained["W1"], trained["b1"]
    W2, b2 = trained["W2"], trained["b2"]
    W3, b3 = trained["W3"], trained["b3"]
    h1_dim = W1.shape[0]
    h2_dim = W2.shape[0]
    n_states = W3.shape[0]

    with open(rust_file, "w") as f:
        f.write("// Auto-generated encoder weights for structural alphabet\n")
        f.write(f"// Trained on {len(pdbs) - n_skipped} structures, "
                f"{len(all_features)} residues\n")
        f.write(f"// Accuracy: {accuracy:.1%} (train), {val_accuracy:.1%} (val)\n")
        f.write(f"// Architecture: 10 → {h1_dim} (ReLU) → {h2_dim} (ReLU) → {n_states}\n")
        f.write(f"// Standardization folded into W1/b1\n\n")

        f.write(f"const H1_DIM: usize = {h1_dim};\n")
        f.write(f"const H2_DIM: usize = {h2_dim};\n\n")

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

        centroids = trained["centroids"]
        ed = trained["embed_dim"]
        f.write(f"static CENTROIDS: [[f64; {ed}]; NUM_STATES] = [\n")
        for row in centroids:
            vals = ", ".join(f"{v: .10f}" for v in row)
            f.write(f"    [{vals}],\n")
        f.write("];\n")

    print(f"  Rust constants saved to {rust_file}")


if __name__ == "__main__":
    main()
