# Ferritin — Next Steps

## Current State (2026-04-02)

- **248 tests** (116 Rust + 132 Python), all passing, zero warnings
- Full PDB/mmCIF I/O with hierarchy navigation + numpy arrays
- 4 alignment algorithms (TM, SOI, Flex, MM) each with single/1-to-N/N-to-M + rayon parallelism
- 6 geometry building blocks (Kabsch, RMSD, apply_transform, secondary structure, TM-score)
- 3-layer architecture: pure Rust → PyO3 connector → pure Python

---

## Priority 1: Oracle Testing

Validate our results against independent implementations. No function untested policy.

- [x] Set up **Biopython** + **Gemmi** as Python-side oracles (pip install)
- [x] Create `tests/oracle/` with shared PDB fixtures (1ubq PDB+CIF, 4hhb, 1crn)
- [x] Oracle tests for I/O: atom counts, chain counts, coordinates (3dp), B-factors, occupancies, elements, atom names, residue names, chain IDs — 72 tests, all matching across 3 tools
- [x] Compare TM-scores against C++ USAlign binary — 7 protein pairs, TM-scores match to 5e-4, RMSD to 0.05 — 35 tests
- [ ] Set up **BALL** as C++ oracle (compiled, Python bindings broken — use via subprocess)
- [ ] Oracle tests for edge cases: multi-model NMR, insertion codes, alt conformers (tools disagree on alt conformer handling — needs specialized tests)
- [ ] CI job that runs oracle tests on a schedule (not every commit — they're slow)

## Priority 2: Selection / Filtering

BALL-inspired string-based selection language for interactive use.

- [ ] Design selection grammar: `"name CA and chain A"`, `"backbone and not hetero"`, `"resid 1:50"`
- [ ] Implement parser in Python layer (not Rust — it's syntactic sugar over filter closures)
- [ ] `structure.select("name CA and chain A")` → filtered `Structure` or atom list
- [ ] Common predicates: `name`, `element`, `chain`, `resname`, `resid`, `backbone`, `sidechain`, `hetero`
- [ ] Boolean combinators: `and`, `or`, `not`, parentheses
- [ ] Consider whether to also expose as Rust API or keep Python-only

## Priority 3: DataFrame Integration

Make structures play nicely with pandas / polars.

- [ ] `structure.to_dataframe()` → pandas DataFrame with columns: atom_name, element, x, y, z, b_factor, occupancy, residue_name, residue_serial, chain_id, hetero
- [ ] `Structure.from_dataframe(df)` → round-trip back (nice to have, lower priority)
- [ ] Residue-level DataFrame: one row per residue with phi/psi angles, SS assignment
- [ ] Alignment result to DataFrame: per-position distances, aligned residue pairs

## Priority 4: Wire CLI mm=1

The `usalign` binary's `-mm 1` mode currently prints "not yet wired".

- [ ] Add `run_mmalign_mode()` in `usalign.rs` that calls `mmalign_complex()`
- [ ] Parse chain lists from PDB, build `ChainData` per chain
- [ ] Output formatting: chain assignments, per-chain TM-scores, alignment strings
- [ ] `-m matrix.txt` rotation matrix output for MM-align
- [ ] Batch mode: `-dir` for all-vs-all complex comparison

## Priority 5: MM-align Refinements

The orchestrator works but has room for improvement.

- [ ] Extract actual transforms from SE results (currently using identity placeholder)
- [ ] Wire dimer detection + `adjust_dimer_assignment()` for 2+2 chain cases
- [ ] Use `homo_refined_greedy_search` with real per-pair rotation matrices (not identity)
- [ ] Use `hetero_refined_greedy_search` when appropriate
- [ ] Add `trim_complex()` preprocessing for very large chains
- [ ] Validate against C++ USAlign `-mm 1` output on reference complexes

## Priority 6: More Building Blocks

Functions useful for the broader structural bioinformatics community.

- [ ] **Superpose and return transformed coords**: `superpose(mobile, reference)` → transformed coords
- [ ] **Contact map**: `contact_map(structure, cutoff=8.0)` → NxN boolean or distance matrix
- [ ] **Phi/Psi angles**: `ramachandran_angles(structure)` → per-residue backbone torsions
- [ ] **SASA**: solvent-accessible surface area (if pdbtbx supports it or implement)
- [ ] **Distance matrix**: `distance_matrix(coords)` → NxN numpy array
- [ ] **Batch RMSD**: `rmsd_matrix(structures, n_threads)` → NxN TM-score or RMSD matrix

## Priority 7: Packaging & Distribution

- [ ] Publish `ferritin-connector` wheel to PyPI via maturin
- [ ] Publish `ferritin` Python package to PyPI
- [ ] CI/CD pipeline (GitHub Actions): build wheels for Linux/macOS/Windows
- [ ] Sphinx or mkdocs documentation site
- [ ] Example notebooks (Jupyter)
- [ ] Benchmark suite vs Biopython/Gemmi/MDAnalysis for I/O and alignment speed

## Ideas (Lower Priority)

- [ ] **CP-align** (circular permutation) wired to Python — already implemented in Rust
- [ ] **NW-align** (sequence alignment) exposed as building block
- [ ] **Batch loading**: `load_many(paths, n_threads)` for parallel PDB loading
- [ ] **mmCIF writing** from modified structures (edit atoms, save back)
- [ ] **R*-tree spatial queries** exposed to Python (pdbtbx has this with rstar feature)
- [ ] **Visualization helpers**: `structure.to_nglview()`, `structure.to_py3dmol()`
- [ ] **BiochemicalAlgorithms.jl-style column access**: `structure.atoms.positions` returning array view
