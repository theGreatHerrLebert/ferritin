"""Fold a structure using proteon features + OpenFold model.

Proves that proteon's supervision/MSA features are consumable by
OpenFold's forward pass. Loads a PDB, extracts features via proteon,
maps them to OpenFold's feature dict format, runs inference with
pre-trained weights, and writes the predicted structure.

Usage:
    python validation/fold_with_proteon_features.py \
        --pdb test-pdbs/1crn.pdb \
        --weights /path/to/finetuning_no_templ_ptm_2.pt \
        --out-pdb predicted_1crn.pdb
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import torch

OPENFOLD_ROOT = Path("/globalscratch/dateschn/proteon-benchmark/openfold")
sys.path.insert(0, str(OPENFOLD_ROOT))

# ---------------------------------------------------------------------------
# AA ordering maps (proteon → OpenFold)
# ---------------------------------------------------------------------------

_OF_RESTYPES = "ARNDCQEGHILKMFPSTWYV"
_OF_AA_TO_IDX = {aa: i for i, aa in enumerate(_OF_RESTYPES)}
_PT_RESTYPES = "ACDEFGHIKLMNPQRSTVWY"
_PT_TO_OF = np.array([_OF_AA_TO_IDX[aa] for aa in _PT_RESTYPES] + [20], dtype=np.int32)


def proteon_aatype_to_openfold(aatype: np.ndarray) -> np.ndarray:
    return _PT_TO_OF[aatype]


# ---------------------------------------------------------------------------
# Feature adapter: proteon → OpenFold feature dict
# ---------------------------------------------------------------------------


def proteon_to_openfold_features(
    structure_ex,
    sequence_ex,
    *,
    num_recycling: int = 3,
    device: torch.device = torch.device("cpu"),
) -> Dict[str, torch.Tensor]:
    """Convert proteon StructureSupervisionExample + SequenceExample
    to the feature dict that OpenFold's AlphaFold.forward() expects.
    """
    n_res = structure_ex.length

    # --- aatype (mapped to OpenFold ordering) ---
    aatype = proteon_aatype_to_openfold(structure_ex.aatype)
    aatype_t = torch.tensor(aatype, dtype=torch.int64, device=device)

    # --- residue_index ---
    residue_index = torch.arange(n_res, dtype=torch.int64, device=device)

    # --- target_feat: [N, 22] = [between_segment(1), aatype_one_hot(21)] ---
    between_segment = torch.zeros(n_res, 1, dtype=torch.float32, device=device)
    aatype_one_hot = torch.nn.functional.one_hot(aatype_t, 21).float()
    target_feat = torch.cat([between_segment, aatype_one_hot], dim=-1)

    # --- MSA features ---
    if sequence_ex is not None and sequence_ex.msa is not None:
        msa_raw = sequence_ex.msa  # (N_seq, N_res) int32, proteon encoding
        n_seq = msa_raw.shape[0]
        # Map MSA from proteon AA order to OpenFold order
        # proteon gap=21, OpenFold gap=21, X=20 in both
        msa_mapped = np.zeros_like(msa_raw)
        for i in range(n_seq):
            for j in range(n_res):
                v = msa_raw[i, j]
                if v < 20:
                    msa_mapped[i, j] = _PT_TO_OF[v]
                else:
                    msa_mapped[i, j] = v  # gap (21) and X (20) stay the same

        msa_t = torch.tensor(msa_mapped, dtype=torch.int64, device=device)
        del_matrix = torch.tensor(
            sequence_ex.deletion_matrix, dtype=torch.float32, device=device
        )
    else:
        # Single-sequence MSA (query only)
        n_seq = 1
        msa_t = aatype_t.unsqueeze(0)  # (1, N_res)
        del_matrix = torch.zeros(1, n_res, dtype=torch.float32, device=device)

    # msa_feat: [N_seq, N_res, 49]
    # = [msa_one_hot(23), has_deletion(1), deletion_value(1),
    #    cluster_profile(23), deletion_mean(1)]
    msa_one_hot = torch.nn.functional.one_hot(msa_t, 23).float()  # (S, N, 23)
    has_deletion = (del_matrix > 0).float().unsqueeze(-1)  # (S, N, 1)
    deletion_value = (
        torch.atan(del_matrix / 3.0) * (2.0 / np.pi)
    ).unsqueeze(-1)  # (S, N, 1)

    # cluster_profile: average AA distribution across MSA
    cluster_profile = msa_one_hot.mean(dim=0, keepdim=True).expand(
        n_seq, -1, -1
    )  # (S, N, 23)
    deletion_mean = del_matrix.mean(dim=0, keepdim=True)
    deletion_mean = (
        torch.atan(deletion_mean / 3.0) * (2.0 / np.pi)
    ).unsqueeze(-1).expand(n_seq, -1, -1)  # (S, N, 1)

    msa_feat = torch.cat(
        [msa_one_hot, has_deletion, deletion_value, cluster_profile, deletion_mean],
        dim=-1,
    )  # (S, N, 49)

    # --- Masks ---
    seq_mask = torch.ones(n_res, dtype=torch.float32, device=device)
    msa_mask = torch.ones(n_seq, n_res, dtype=torch.float32, device=device)

    # --- atom37_atom_exists ---
    atom37_atom_exists = torch.tensor(
        structure_ex.atom37_atom_exists, dtype=torch.int64, device=device
    )

    # --- Extra MSA (empty — we put everything in the main MSA) ---
    extra_msa = torch.zeros(1, n_res, dtype=torch.int64, device=device)
    extra_has_deletion = torch.zeros(1, n_res, dtype=torch.float32, device=device)
    extra_deletion_value = torch.zeros(1, n_res, dtype=torch.float32, device=device)
    extra_msa_mask = torch.zeros(1, n_res, dtype=torch.float32, device=device)

    # --- Atom index mappings (needed for final atom14→atom37 conversion) ---
    residx_atom14_to_atom37 = torch.tensor(
        structure_ex.residx_atom14_to_atom37, dtype=torch.int64, device=device
    )
    residx_atom37_to_atom14 = torch.tensor(
        structure_ex.residx_atom37_to_atom14, dtype=torch.int64, device=device
    )
    atom14_atom_exists = torch.tensor(
        structure_ex.atom14_atom_exists, dtype=torch.float32, device=device
    )

    # --- Ground truth for loss/metrics ---
    all_atom_positions = torch.tensor(
        structure_ex.all_atom_positions, dtype=torch.float32, device=device
    )
    all_atom_mask = torch.tensor(
        structure_ex.all_atom_mask, dtype=torch.float32, device=device
    )

    # --- Assemble feature dict ---
    feats = {
        "aatype": aatype_t,
        "target_feat": target_feat,
        "residue_index": residue_index,
        "msa_feat": msa_feat,
        "seq_mask": seq_mask,
        "msa_mask": msa_mask,
        "atom37_atom_exists": atom37_atom_exists,
        "extra_msa": extra_msa,
        "extra_has_deletion": extra_has_deletion,
        "extra_deletion_value": extra_deletion_value,
        "extra_msa_mask": extra_msa_mask,
        "residx_atom14_to_atom37": residx_atom14_to_atom37,
        "residx_atom37_to_atom14": residx_atom37_to_atom14,
        "atom14_atom_exists": atom14_atom_exists,
        "all_atom_positions": all_atom_positions,
        "all_atom_mask": all_atom_mask,
        "no_recycling_iters": torch.tensor(num_recycling, dtype=torch.int64, device=device),
    }

    # Add recycling dimension: stack num_recycling+1 copies along last dim.
    # OpenFold's forward() determines num_iters from aatype.shape[-1] and
    # slices every tensor with [..., cycle_no]. no_recycling_iters is not
    # part of the batch dict — it's only used by the feature pipeline.
    recycling_feats = {}
    for k, v in feats.items():
        if k == "no_recycling_iters":
            continue  # not part of the model batch
        recycling_feats[k] = torch.stack([v] * (num_recycling + 1), dim=-1)

    return recycling_feats


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Fold with proteon features + OpenFold")
    parser.add_argument("--pdb", type=Path, required=True, help="Input PDB file")
    parser.add_argument("--weights", type=Path, required=True, help="OpenFold checkpoint")
    parser.add_argument("--config", type=str, default="finetuning_no_templ_ptm",
                        help="OpenFold config preset")
    parser.add_argument("--out-pdb", type=Path, default=None, help="Output PDB path")
    parser.add_argument("--device", type=str, default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--msa-db", type=Path, default=None, help="MSA DB prefix for search")
    parser.add_argument("--msa-kmi", type=Path, default=None, help="KMI path for MSA search")
    args = parser.parse_args()

    device = torch.device(args.device)
    print(f"Device: {device}")
    torch.set_grad_enabled(False)
    if torch.cuda.is_available():
        torch.set_float32_matmul_precision("high")

    # --- Load structure with proteon ---
    print(f"\nLoading {args.pdb} with proteon ...")
    import proteon

    structure = proteon.load(str(args.pdb))
    chain = structure.chains[0]
    print(f"  Chain {chain.id}, {len([r for r in chain.residues if r.is_amino_acid])} residues")

    structure_ex = proteon.build_structure_supervision_example(structure, chain_id=chain.id)
    print(f"  Structure supervision: {structure_ex.length} residues, all fields present")

    # --- Optional MSA ---
    sequence_ex = None
    if args.msa_db and args.msa_kmi:
        print(f"\nSearching MSA against {args.msa_db} ...")
        from proteon.msa_backend import (
            open_search_engine_from_mmseqs_db_with_kmi,
            build_sequence_example_with_msa,
        )
        engine = open_search_engine_from_mmseqs_db_with_kmi(
            str(args.msa_db), str(args.msa_kmi),
        )
        sequence_ex = build_sequence_example_with_msa(
            structure, engine, chain_id=chain.id, max_seqs=256,
        )
        print(f"  MSA depth: {sequence_ex.msa.shape[0]}")
    else:
        print("\n  No MSA DB provided — running single-sequence mode")

    # --- Build OpenFold features ---
    print("\nBuilding OpenFold feature dict ...")
    feats = proteon_to_openfold_features(structure_ex, sequence_ex, device=device)
    for k, v in feats.items():
        if isinstance(v, torch.Tensor):
            print(f"  {k}: {v.shape} {v.dtype}")

    # --- Load model ---
    print(f"\nLoading OpenFold model ({args.config}) ...")
    from openfold.config import model_config
    from openfold.model.model import AlphaFold
    from openfold.utils.import_weights import import_openfold_weights_

    config = model_config(args.config)
    model = AlphaFold(config)
    model = model.eval()

    checkpoint = torch.load(str(args.weights), map_location="cpu", weights_only=False)
    if "ema" in checkpoint:
        state_dict = checkpoint["ema"]["params"]
    elif "model_state_dict" in checkpoint:
        state_dict = checkpoint["model_state_dict"]
    else:
        state_dict = checkpoint

    import_openfold_weights_(model=model, state_dict=state_dict)
    model = model.to(device)
    print(f"  Model loaded, {sum(p.numel() for p in model.parameters()) / 1e6:.1f}M parameters")

    # --- Forward pass ---
    print("\nRunning forward pass ...")
    t0 = time.monotonic()
    with torch.no_grad():
        outputs = model(feats)
    dt = time.monotonic() - t0
    print(f"  Forward pass: {dt:.1f}s")

    # --- Extract predicted structure ---
    final_positions = outputs["final_atom_positions"]  # (N, 37, 3)
    final_mask = outputs["final_atom_mask"]  # (N, 37)
    print(f"  Predicted positions: {final_positions.shape}")

    # pLDDT confidence
    if "plddt" in outputs:
        plddt = outputs["plddt"].cpu().numpy()  # (N,) or (N, 37)
        if plddt.ndim == 2:
            plddt_per_res = plddt[:, 1]  # CA
        else:
            plddt_per_res = plddt
        print(f"  pLDDT: mean={plddt_per_res.mean():.1f}, "
              f"median={np.median(plddt_per_res):.1f}, "
              f"min={plddt_per_res.min():.1f}")

    # ptm score
    if "predicted_tm_score" in outputs:
        ptm = outputs["predicted_tm_score"].item()
        print(f"  pTM: {ptm:.3f}")

    # --- Write output PDB ---
    if args.out_pdb:
        from openfold.np import protein as protein_module

        positions_np = final_positions.cpu().numpy()
        mask_np = final_mask.cpu().numpy()
        aatype_np = proteon_aatype_to_openfold(structure_ex.aatype)

        # pLDDT as b-factors
        if "plddt" in outputs:
            b_factors = np.zeros_like(mask_np)
            plddt_np = outputs["plddt"].cpu().numpy()
            if plddt_np.ndim == 2:
                b_factors = plddt_np
            else:
                b_factors = np.repeat(plddt_np[:, None], 37, axis=1)
        else:
            b_factors = np.zeros_like(mask_np)

        prot = protein_module.Protein(
            atom_positions=positions_np,
            aatype=aatype_np,
            atom_mask=mask_np,
            residue_index=np.arange(structure_ex.length),
            b_factors=b_factors,
            chain_index=np.zeros(structure_ex.length, dtype=np.int32),
        )
        pdb_str = protein_module.to_pdb(prot)
        args.out_pdb.write_text(pdb_str)
        print(f"\n  Wrote predicted structure to {args.out_pdb}")

    print("\nDone!")


if __name__ == "__main__":
    main()
