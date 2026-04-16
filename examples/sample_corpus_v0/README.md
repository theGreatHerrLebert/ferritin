# sample_corpus_v0

A reference artifact tree produced by `build_local_corpus_smoke_release` on 10
PDBs from `/scratch/TMAlign/test-pdbs/` (33 chains after expansion).

Built 2026-04-16 with `split_ratios={train: 0.7, val: 0.15, test: 0.15}`.

## Layout

```
sample_corpus_v0/
├── corpus/
│   ├── corpus_release_manifest.json   # top-level manifest (format: ferritin.corpus_release.v0)
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
    ├── release_manifest.json          # includes split_counts
    └── training_examples.jsonl        # one row per chain-record with `split` field
```

## What's missing from the committed tree

`tensors.npz` files are stripped (they dominate size at ~3.7 MB each).
To regenerate them (and the whole tree) with the same inputs:

```bash
python examples/10_corpus_release_smoke.py \
    --out sample_corpus_v0_full \
    /scratch/TMAlign/test-pdbs/1a*.pdb /scratch/TMAlign/test-pdbs/1b*.pdb
```

Manifest `tensor_sha256` entries are retained as a record of what was produced.
They refer to the stripped `tensors.npz`; re-running locally may or may not
reproduce the exact hash depending on floating-point ordering of the build.
