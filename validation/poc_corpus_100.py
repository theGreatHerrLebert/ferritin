"""Proof-of-concept: full corpus release with MSA on 100 PDBs.

Runs the complete AF2-style training data pipeline:
  load → prep → supervise → MSA search → sequence features → join → release

Validates feature completeness and correctness at small scale before
committing to archive-scale runs.

Usage:
    source venv/bin/activate
    python validation/poc_corpus_100.py
"""

from __future__ import annotations

import random
import sys
import time
from pathlib import Path

PDB_DIR = Path("/globalscratch/dateschn/proteon-benchmark/pdbs_50k")
DB_PREFIX = Path("/globalscratch/dateschn/uniref_subset/uniref50_1m")
KMI_PATH = Path("/globalscratch/dateschn/uniref_subset/uniref50_1m_ext.kmi")
OUT_DIR = Path("/globalscratch/dateschn/proteon-benchmark/poc_corpus_100")
N_STRUCTURES = 100
SEED = 42
CHUNK_SIZE = 25  # process 25 structures at a time to bound memory


def main():
    t_start = time.monotonic()

    # --- Sample 100 PDBs ---
    all_pdbs = sorted(PDB_DIR.glob("*.pdb"))
    print(f"Available PDBs: {len(all_pdbs)}")

    rng = random.Random(SEED)
    sample = rng.sample(all_pdbs, min(N_STRUCTURES, len(all_pdbs)))
    sample.sort()
    print(f"Sampled: {len(sample)} PDBs")
    print(f"First 5: {[p.stem for p in sample[:5]]}")

    # --- Open MSA engine ---
    print(f"\nOpening MSA engine from {DB_PREFIX} + {KMI_PATH} ...")
    from proteon.msa_backend import open_search_engine_from_mmseqs_db_with_kmi

    t0 = time.monotonic()
    engine = open_search_engine_from_mmseqs_db_with_kmi(
        str(DB_PREFIX), str(KMI_PATH),
    )
    print(f"Engine ready ({time.monotonic() - t0:.1f}s)")

    # --- Run the full corpus release pipeline ---
    print(f"\nBuilding corpus release at {OUT_DIR} ...")
    print(f"  chunk_size={CHUNK_SIZE}")
    print(f"  msa_max_seqs=256")
    print()

    from proteon.corpus_smoke import build_local_corpus_smoke_release

    t0 = time.monotonic()
    result = build_local_corpus_smoke_release(
        paths=sample,
        out_dir=OUT_DIR,
        release_id="poc-100-msa",
        code_rev=None,
        config_rev=None,
        msa_engine=engine,
        msa_max_seqs=256,
        chunk_size=CHUNK_SIZE,
        rescue_load=True,
        overwrite=True,
    )
    pipeline_time = time.monotonic() - t0

    # --- Report ---
    import json
    import numpy as np

    print(f"\n{'='*60}")
    print(f"PoC Corpus Release Complete")
    print(f"{'='*60}")
    print(f"Output: {result}")
    print(f"Pipeline time: {pipeline_time:.1f}s ({pipeline_time/60:.1f} min)")

    # Read the corpus manifest
    manifest_path = OUT_DIR / "corpus" / "corpus_release_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        print(f"\nCorpus manifest:")
        print(f"  Release ID: {manifest.get('release_id')}")
        counts = manifest.get("counts", {})
        for k, v in counts.items():
            print(f"  {k}: {v}")

    # Read the validation report
    val_path = OUT_DIR / "corpus" / "validation_report.json"
    if val_path.exists():
        val = json.loads(val_path.read_text())
        print(f"\nValidation: {'PASS' if val.get('valid') else 'FAIL'}")
        if not val.get("valid"):
            for err in val.get("errors", [])[:5]:
                print(f"  ERROR: {err}")

    # Check supervision tensors
    sup_manifest_path = OUT_DIR / "prepared" / "supervision_release" / "release_manifest.json"
    if sup_manifest_path.exists():
        sup = json.loads(sup_manifest_path.read_text())
        print(f"\nStructure supervision:")
        print(f"  Examples: {sup.get('count_examples', 0)}")
        print(f"  Failures: {sup.get('count_failures', 0)}")

    # Check sequence tensors
    seq_manifest_path = OUT_DIR / "sequence" / "release_manifest.json"
    if seq_manifest_path.exists():
        seq_m = json.loads(seq_manifest_path.read_text())
        print(f"\nSequence features:")
        print(f"  Examples: {seq_m.get('count_examples', 0)}")
        print(f"  Failures: {seq_m.get('count_failures', 0)}")

    # Spot-check: load a few training examples and verify field presence
    print(f"\n--- Spot-check: loading training examples ---")
    try:
        from proteon.supervision_export import load_structure_supervision_examples
        from proteon.sequence_export import load_sequence_examples

        sup_examples = load_structure_supervision_examples(
            OUT_DIR / "prepared" / "supervision_release" / "examples",
            verify_checksum=False,
        )
        print(f"Loaded {len(sup_examples)} structure supervision examples")
        if sup_examples:
            ex = sup_examples[0]
            print(f"\n  Example: {ex.record_id} (chain={ex.chain_id}, length={ex.length})")
            print(f"  Fields present:")
            for field_name in [
                "aatype", "all_atom_positions", "all_atom_mask",
                "atom14_gt_positions", "atom14_alt_gt_positions", "atom14_alt_gt_exists",
                "atom14_atom_is_ambiguous",
                "pseudo_beta", "chi_angles", "rigidgroups_gt_frames",
            ]:
                val = getattr(ex, field_name, None)
                if val is not None:
                    print(f"    {field_name}: {val.shape} {val.dtype}")
                else:
                    print(f"    {field_name}: MISSING")

        seq_examples = load_sequence_examples(
            OUT_DIR / "sequence" / "examples",
            verify_checksum=False,
        )
        print(f"\nLoaded {len(seq_examples)} sequence examples")
        if seq_examples:
            se = seq_examples[0]
            print(f"  Example: {se.record_id} (length={se.length})")
            has_msa = se.msa is not None
            print(f"  MSA present: {has_msa}")
            if has_msa:
                print(f"    msa shape: {se.msa.shape}")
                print(f"    deletion_matrix shape: {se.deletion_matrix.shape}")
                print(f"    msa_mask shape: {se.msa_mask.shape}")
            has_profile = se.msa_profile is not None
            print(f"  MSA profile present: {has_profile}")
            if has_profile:
                print(f"    msa_profile shape: {se.msa_profile.shape}")

            # Count how many sequence examples have MSA
            n_with_msa = sum(1 for s in seq_examples if s.msa is not None and s.msa.shape[0] > 1)
            print(f"\n  Sequences with MSA (depth > 1): {n_with_msa}/{len(seq_examples)}")

    except Exception as e:
        print(f"  Spot-check failed: {e}")

    total_time = time.monotonic() - t_start
    print(f"\nTotal wall time: {total_time:.1f}s ({total_time/60:.1f} min)")


if __name__ == "__main__":
    main()
