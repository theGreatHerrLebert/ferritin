"""Batch fold: proteon features → OpenFold → predicted structures.

Folds N structures from a PDB directory, optionally with MSA, and
reports pLDDT + TM-score against the crystal structure.

Usage:
    python validation/batch_fold.py \
        --pdb-dir /globalscratch/dateschn/proteon-benchmark/pdbs_50k \
        --n 250 --seed 42 \
        --msa-db /globalscratch/dateschn/uniref_subset/uniref50_1m \
        --msa-kmi /globalscratch/dateschn/uniref_subset/uniref50_1m_ext.kmi \
        --out-dir /globalscratch/dateschn/proteon-benchmark/batch_fold_250
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import time
from pathlib import Path

import numpy as np
import torch

OPENFOLD_ROOT = Path("/globalscratch/dateschn/proteon-benchmark/openfold")
PROTEON_VALIDATION = Path(__file__).resolve().parent
sys.path.insert(0, str(OPENFOLD_ROOT))
sys.path.insert(0, str(PROTEON_VALIDATION.parent))

from validation.fold_with_proteon_features import (
    proteon_to_openfold_features,
    proteon_aatype_to_openfold,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-dir", type=Path, required=True)
    parser.add_argument("--n", type=int, default=250)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--weights", type=Path,
                        default=Path("/globalscratch/dateschn/proteon-benchmark/openfold_weights/models--nz--OpenFold/snapshots/6344b2a680cbbae20b8afdb7f154228cc5b0cb27/finetuning_no_templ_ptm_1.pt"))
    parser.add_argument("--config", type=str, default="finetuning_no_templ_ptm")
    parser.add_argument("--msa-db", type=Path, default=None)
    parser.add_argument("--msa-kmi", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--device", type=str, default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--max-length", type=int, default=400,
                        help="Skip chains longer than this (GPU memory)")
    args = parser.parse_args()

    device = torch.device(args.device)
    torch.set_grad_enabled(False)
    if torch.cuda.is_available():
        torch.set_float32_matmul_precision("high")

    # Sample PDBs
    all_pdbs = sorted(args.pdb_dir.glob("*.pdb"))
    rng = random.Random(args.seed)
    sample = rng.sample(all_pdbs, min(args.n, len(all_pdbs)))
    sample.sort()
    print(f"Sampled {len(sample)} PDBs")

    # Load model
    print(f"Loading model ({args.config}) ...")
    from openfold.config import model_config
    from openfold.model.model import AlphaFold
    from openfold.utils.import_weights import import_openfold_weights_

    config = model_config(args.config)
    model = AlphaFold(config).eval()
    ckpt = torch.load(str(args.weights), map_location="cpu", weights_only=False)
    state_dict = ckpt.get("ema", {}).get("params", ckpt)
    import_openfold_weights_(model=model, state_dict=state_dict)
    model = model.to(device)
    print(f"Model ready ({sum(p.numel() for p in model.parameters()) / 1e6:.1f}M params)")

    # MSA engine
    msa_engine = None
    if args.msa_db and args.msa_kmi:
        print(f"Opening MSA engine ...")
        from proteon.msa_backend import open_search_engine_from_mmseqs_db_with_kmi
        msa_engine = open_search_engine_from_mmseqs_db_with_kmi(
            str(args.msa_db), str(args.msa_kmi))
        print("MSA engine ready")

    # Output
    if args.out_dir:
        args.out_dir.mkdir(parents=True, exist_ok=True)

    import proteon
    from openfold.np import protein as protein_module

    results = []
    t_start = time.monotonic()

    for i, pdb_path in enumerate(sample):
        t0 = time.monotonic()
        try:
            structure = proteon.load(str(pdb_path))
            chain = structure.chains[0]
            residues = [r for r in chain.residues if r.is_amino_acid]
            n_res = len(residues)

            if n_res < 10:
                print(f"  [{i+1}/{len(sample)}] {pdb_path.stem}: skip (too short, {n_res} res)")
                continue
            if n_res > args.max_length:
                print(f"  [{i+1}/{len(sample)}] {pdb_path.stem}: skip (too long, {n_res} > {args.max_length})")
                continue

            structure_ex = proteon.build_structure_supervision_example(
                structure, chain_id=chain.id)

            # MSA
            sequence_ex = None
            if msa_engine is not None:
                from proteon.msa_backend import build_sequence_example_with_msa
                sequence_ex = build_sequence_example_with_msa(
                    structure, msa_engine, chain_id=chain.id, max_seqs=256)

            # Build features & fold
            feats = proteon_to_openfold_features(
                structure_ex, sequence_ex, device=device)
            with torch.no_grad():
                outputs = model(feats)

            # Extract results
            positions = outputs["final_atom_positions"].cpu().numpy()
            mask = outputs["final_atom_mask"].cpu().numpy()
            plddt = outputs["plddt"].cpu().numpy()
            if plddt.ndim == 2:
                plddt_mean = float(plddt[:, 1].mean())  # CA
            else:
                plddt_mean = float(plddt.mean())

            ptm = outputs.get("predicted_tm_score")
            ptm_val = float(ptm.item()) if ptm is not None else None

            # Write predicted PDB
            if args.out_dir:
                aatype_of = proteon_aatype_to_openfold(structure_ex.aatype)
                b_factors = plddt if plddt.ndim == 2 else np.repeat(plddt[:, None], 37, axis=1)
                prot = protein_module.Protein(
                    atom_positions=positions,
                    aatype=aatype_of,
                    atom_mask=mask,
                    residue_index=np.arange(n_res),
                    b_factors=b_factors,
                    chain_index=np.zeros(n_res, dtype=np.int32),
                )
                out_pdb = args.out_dir / f"{pdb_path.stem}_pred.pdb"
                out_pdb.write_text(protein_module.to_pdb(prot))

            # TM-score against crystal
            tm1 = None
            try:
                if args.out_dir:
                    pred = proteon.load(str(out_pdb))
                    r = proteon.tm_align(structure, pred)
                    tm1 = float(r.tm_score_chain1)
            except Exception:
                pass

            dt = time.monotonic() - t0
            msa_depth = sequence_ex.msa.shape[0] if sequence_ex and sequence_ex.msa is not None else 0

            result = {
                "pdb": pdb_path.stem,
                "length": n_res,
                "msa_depth": msa_depth,
                "plddt": round(plddt_mean, 1),
                "ptm": round(ptm_val, 3) if ptm_val else None,
                "tm_score": round(tm1, 4) if tm1 else None,
                "time_s": round(dt, 1),
            }
            results.append(result)

            tm_str = f"TM={tm1:.3f}" if tm1 else "TM=N/A"
            ptm_str = f"pTM={ptm_val:.3f}" if ptm_val is not None else "pTM=N/A"
            print(f"  [{i+1}/{len(sample)}] {pdb_path.stem}: L={n_res} pLDDT={plddt_mean:.1f} "
                  f"{ptm_str} {tm_str} MSA={msa_depth} ({dt:.1f}s)")

        except Exception as e:
            dt = time.monotonic() - t0
            print(f"  [{i+1}/{len(sample)}] {pdb_path.stem}: ERROR ({e}) ({dt:.1f}s)")
            results.append({"pdb": pdb_path.stem, "error": str(e)})

    # Summary
    total_time = time.monotonic() - t_start
    successful = [r for r in results if "plddt" in r]
    tm_scores = [r["tm_score"] for r in successful if r.get("tm_score") is not None]
    plddts = [r["plddt"] for r in successful]

    print(f"\n{'='*60}")
    print(f"Batch Fold Summary")
    print(f"{'='*60}")
    print(f"Total: {len(sample)}  Folded: {len(successful)}  Errors: {len(results) - len(successful)}")
    print(f"Wall time: {total_time:.0f}s ({total_time/60:.1f} min)")
    if plddts:
        arr = np.array(plddts)
        print(f"pLDDT:    mean={arr.mean():.1f}  median={np.median(arr):.1f}  "
              f"min={arr.min():.1f}  max={arr.max():.1f}")
    if tm_scores:
        arr = np.array(tm_scores)
        print(f"TM-score: mean={arr.mean():.3f}  median={np.median(arr):.3f}  "
              f"min={arr.min():.3f}  max={arr.max():.3f}")
        print(f"TM>0.5:   {(arr > 0.5).sum()}/{len(arr)} ({100*(arr > 0.5).mean():.1f}%)")
        print(f"TM>0.7:   {(arr > 0.7).sum()}/{len(arr)} ({100*(arr > 0.7).mean():.1f}%)")
        print(f"TM>0.9:   {(arr > 0.9).sum()}/{len(arr)} ({100*(arr > 0.9).mean():.1f}%)")

    if args.out_dir:
        report_path = args.out_dir / "fold_report.json"
        report = {
            "config": {
                "n": len(sample), "seed": args.seed,
                "msa": args.msa_db is not None,
                "max_length": args.max_length,
                "device": str(device),
            },
            "summary": {
                "n_folded": len(successful),
                "n_errors": len(results) - len(successful),
                "wall_time_s": round(total_time, 1),
                "mean_plddt": round(float(np.mean(plddts)), 1) if plddts else None,
                "median_tm": round(float(np.median(tm_scores)), 3) if tm_scores else None,
            },
            "per_structure": results,
        }
        report_path.write_text(json.dumps(report, indent=2))
        print(f"\nReport: {report_path}")


if __name__ == "__main__":
    main()
