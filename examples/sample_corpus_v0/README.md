# sample_corpus_v0

A reference artifact tree produced by `build_local_corpus_smoke_release` on 10
PDBs from `/scratch/TMAlign/test-pdbs/` (33 chain records after multi-chain
expansion). Zero failures, all three splits populated.

Built 2026-04-16 with `split_ratios={train: 0.7, val: 0.15, test: 0.15}`.

## Layout

```
sample_corpus_v0/
├── corpus/
│   ├── corpus_release_manifest.json   # top-level (format: ferritin.corpus_release.v0)
│   └── validation_report.json         # post-build validation (ok=true)
├── prepared/
│   ├── prepared_structures.jsonl      # one row per prepared input
│   └── supervision_release/
│       ├── release_manifest.json
│       ├── failures.jsonl
│       └── examples/
│           ├── manifest.json          # AF2-contract tensor inventory + sha256
│           └── examples.jsonl         # one row per chain-record
├── sequence/
│   ├── release_manifest.json
│   ├── failures.jsonl
│   └── examples/
│       ├── manifest.json
│       └── examples.jsonl
└── training/
    ├── release_manifest.json          # format: ferritin.training_example.parquet.v0
    └── training_examples.jsonl        # one row per chain-record with `split` field
```

## What's missing from the committed tree

Bulk tensor files are stripped so the committed tree stays small and GitHub-
browseable:

- `prepared/supervision_release/examples/tensors.npz` — 3.7 MB
- `sequence/examples/tensors.npz` — 8 KB
- `training/training.parquet` — 5.1 MB

Manifest `*_sha256` entries are retained as a record of what was produced.

To regenerate the full tree with the same inputs:

```bash
python examples/10_corpus_release_smoke.py \
    --out sample_corpus_v0_full \
    /scratch/TMAlign/test-pdbs/1a*.pdb /scratch/TMAlign/test-pdbs/1b*.pdb
```

Re-running may not reproduce the exact hashes depending on floating-point
ordering of the build and Parquet writer-side compression details. The
commits should still pass the release validator because it checks
completeness and per-row ordering, not byte equality to the committed
manifests.
