# CLAUDE.md — Ferritin

Structural bioinformatics toolkit in Rust with Python bindings.

## Repository Overview

### Workspace Crates

| Crate | Purpose |
|-------|---------|
| `ferritin-align/` | Core alignment algorithms (TM-align + US-align), library only |
| `ferritin-io/` | Structure I/O for PDB/mmCIF via pdbtbx |
| `ferritin-bin/` | CLI binaries (`tmalign`, `usalign`) |
| `ferritin-connector/` | PyO3 Python bindings (cdylib) |

### Python Package

| Package | Purpose |
|---------|---------|
| `packages/ferritin/` | Pythonic wrapper over ferritin-connector |

### Reference C++ (read-only)

| Directory | Purpose |
|-----------|---------|
| `TMAlign/` | Original C++ TM-align source |
| `USAlign/` | Original C++ US-align source (multi-file) |

### External Dependencies

| Crate | Path | Purpose |
|-------|------|---------|
| `pdbtbx` | `../pdbtbx` | PDB/mmCIF parsing library |

## Build Commands

```bash
cd ferritin

# Build everything
cargo build

# Build specific crate
cargo build -p ferritin-align

# Run tests
cargo test
cargo test -p ferritin-align    # 95 tests

# Run binaries
cargo run --bin tmalign -- structure1.pdb structure2.pdb
cargo run --bin usalign -- structure1.pdb structure2.pdb
cargo run --bin usalign -- structure1.pdb structure2.pdb --outfmt 2
cargo run --bin usalign -- structure1.pdb structure2.pdb --mm 1   # multi-chain

# Release build
cargo build --release

# PyO3 connector (requires maturin)
cd ferritin-connector
maturin develop --release
```

## Architecture

### ferritin-align Module Layout

```
src/
├── lib.rs           # pub mod core; pub mod ext;
├── core/            # Basic TM-align algorithms
│   ├── types.rs     # Coord3D, Transform, TMParams, AlignResult, AlignOptions, etc.
│   ├── kabsch.rs    # Kabsch optimal rotation (SVD-based)
│   ├── tmscore.rs   # TM-score computation (score_fun8, tmscore8_search, etc.)
│   ├── nwdp.rs      # Needleman-Wunsch DP (4 variants)
│   ├── secondary_structure.rs  # HETC assignment from CA distances
│   ├── residue_map.rs         # AA 3-letter <-> 1-letter conversion
│   └── align/       # Alignment strategies
│       ├── tmalign.rs    # Main TM-align orchestrator
│       ├── cpalign.rs    # Circular permutation
│       ├── dp_iter.rs    # DP iteration refinement
│       └── initial_*.rs  # 5 initialization strategies
└── ext/             # US-align extensions
    ├── blosum.rs    # BLOSUM62 + BLASTN scoring matrix
    ├── nwalign.rs   # Gotoh affine-gap alignment
    ├── se.rs        # Structure Extension refinement
    ├── hwrmsd.rs    # Hybrid weighted RMSD iteration
    ├── flexalign.rs # Flexible hinge-based alignment
    ├── soialign/    # Sequence Order Independent alignment
    │   ├── close_k.rs, sec_bond.rs, greedy.rs, iter.rs
    └── mmalign/     # Multi-chain complex alignment
        ├── chain_assign.rs, complex_score.rs, dimer.rs, iter.rs, trim.rs
```

### PyO3 Pattern (rustims-style)

Three-layer architecture:
1. **Pure Rust** (`ferritin-align`, `ferritin-io`) — no Python dependency
2. **PyO3 connector** (`ferritin-connector`) — `#[pyclass]` wrappers with `inner: RustType`
3. **Python package** (`packages/ferritin/`) — `RustWrapperObject` ABC, Pythonic API

## Key Design Choices

- `core::` = TM-align port, `ext::` = US-align extensions on top
- `DPWorkspace` pre-allocates DP matrices to avoid hot-loop allocations
- `nwdp_core()` is a generic DP engine parameterized by scoring closure
- All functions return values (no mutation of output parameters)
- Float operation order preserved from C++ for numerical fidelity (~4-5 decimal places)
- `d0` formula uses `.powf(1.0/3.0)` (not `.cbrt()`) to match C++ `-ffast-math` behavior

## Numerical Precision

C++ TMalign/USalign are compiled with `-ffast-math`. Rust does not allow this. TM-scores match to ~4-5 decimal places. Key constants: Kabsch `epsilon=1e-8`, `tol=0.01`, `sqrt3=1.73205080756888`.
