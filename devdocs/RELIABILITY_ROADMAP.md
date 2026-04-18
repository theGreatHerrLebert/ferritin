# Ferritin Reliability Roadmap

**Last updated: 2026-04-09**

Ferritin already has a strong compute core, but the reliability story is uneven:

- Core I/O, alignment, Arrow, and Python API tests are in good shape.
- Oracle checks exist and pass locally, but are not part of default CI.
- Structure preparation is partly implemented in Rust internals, but not yet surfaced as a coherent public workflow.
- Validation artifacts exist, but too much of the reliability story still lives in scripts and notes instead of enforced gates.

This document focuses on one goal:

**Make ferritin the default reliable interface for AI-assisted structure-based bioinformatics, especially for data wrangling and structure preparation before downstream modeling.**

---

## Current Baseline

Observed locally on 2026-04-19:

- `cargo test --workspace` passes
- `pytest tests/ -m "not slow"` passes (unit + library-parity oracle
  tests; the library-parity slice of `tests/oracle/` is now exercised
  inline rather than excluded via `-k "not oracle"`)
- `pytest tests/oracle/` passes (library-parity tests + self-skip
  handling when heavier oracles like reduce / USAlign / BALL Julia /
  OpenMM aren't installed)
- `pytest -m slow` covers the 19 >10s tests (minimizer sweeps, corpus
  pipelines); skip in the inner loop with `pytest -m "not slow"`.

Current strengths:

- Pure Rust kernels with thin Python wrappers
- Real oracle tests for I/O and TM-align
- Batch-first API design
- Arrow/Parquet path for scalable feature/data workflows
- ForceField trait with two implementations (AMBER96, CHARMM19+EEF1)
- Full MD pipeline (Verlet, SHAKE, EEF1 solvation, trajectory output)
- Three minimizers (SD, CG, L-BFGS)
- 5K structure validation (run 2: 0 failures, all edge cases triaged)

Current weaknesses:

- CI only runs Ubuntu + Python 3.12
- Default CI excludes oracle tests
- No reliability matrix across OS / Python versions
- No scheduled validation job
- No coverage thresholds
- No end-to-end `prepare()` API for real-world structure cleanup
- No strong regression corpus for ugly archive cases
- CLI/data-engine paths are lightly guarded compared with the Python API
- Force field energy not yet oracle-validated against BALL Julia systematically
- Neighbor list not yet wired into default MD/minimizer paths

---

## Reliability Principles

1. Every public feature should have at least one fast deterministic test.
2. Every scientific claim should have an oracle or invariant behind it.
3. Validation results do not count unless they are reproducible and rerunnable.
4. The default CI path should catch common regressions without relying on local setup.
5. Structure preparation is a product surface, not an internal implementation detail.
6. Data wrangling paths are first-class reliability targets, not just model-prep helpers.

---

## Agent Notes Policy

Ferritin’s boundary-layer docstrings can carry short structured notes for AI-assisted callers, but they should be treated as part of the interface contract, not as free-form commentary.

Allowed prefixes:

- `WATCH`
- `PREFER`
- `COST`
- `INVARIANT`

Rules:

1. Use Agent Notes only on public boundary APIs where misuse is likely.
2. Do not add them to every function by default.
3. Do not repeat obvious signature or return-value information.
4. Prefer notes that affect decision-making:
   scaling, asymmetry, idempotency, interpretation traps, safer batch alternatives, or memory cost.
5. Any note that describes a guarantee or important behavioral caveat should be backed by tests, a stable invariant, or both.
6. README claims about Agent Notes should match the real coverage in the codebase.

---

## Reliability Tiers

### Tier 0: Build and Smoke

Purpose:
- Catch broken imports, packaging regressions, linker issues, and obvious API breaks fast.

Must cover:
- `cargo test --workspace`
- Python package import
- connector build via `maturin develop`
- public API smoke import
- CLI binaries start and show help

Success criteria:
- Runs on every PR
- finishes quickly
- no external oracle dependencies

### Tier 1: Deterministic Functional Tests

Purpose:
- Verify stable behavior of core algorithms and wrapper behavior.

Must cover:
- I/O
- alignment
- analysis
- selection
- SASA
- DSSP
- H-bonds
- forcefield/minimization API contracts
- Arrow/Parquet roundtrip
- CLI integration for `ferritin-ingest`, `tmalign`, `usalign`

Success criteria:
- deterministic
- same assertions on CI and local dev
- no dependence on internet or large downloads

### Tier 2: Oracle Validation

Purpose:
- Compare ferritin against independent trusted tools.

Must cover:
- I/O vs Gemmi / Biopython
- TM-align vs C++ USAlign
- SASA vs FreeSASA
- hydrogen placement / reconstruction vs BALL or BiochemicalAlgorithms.jl
- structure preparation vs PDBFixer where applicable
- DSSP vs mkdssp or GROMACS `gmx dssp`

Success criteria:
- executable in CI on schedule
- thresholds are explicit
- failures are actionable, not advisory

### Tier 3: Corpus and Regression Validation

Purpose:
- Prevent old bugs from reappearing on real archive data.

Must cover:
- multi-model NMR
- altloc-heavy structures
- insertion codes
- missing atoms
- malformed but common PDB quirks
- mixed protein/ligand/water systems
- mmCIF edge cases
- large structures

Success criteria:
- fixed corpus in-repo or pinned externally
- every bug gets a regression fixture
- tolerances are documented

### Tier 4: Performance and Scaling

Purpose:
- Protect ferritin’s value proposition: modern, easy, fast.

Must cover:
- batch loading throughput
- SASA throughput
- alignment throughput
- Arrow ingest throughput
- memory growth on large structures / large corpora

Success criteria:
- benchmark baselines stored
- regressions called out automatically
- no silent algorithmic blowups

---

## Priority Roadmap

## P0: Make Reliability Enforced

These are the highest-priority changes.

### 1. Expand CI Matrix

Add:

- Linux + macOS CI
- Python 3.10, 3.11, 3.12, 3.13 for the Python package
- at least one release-like wheel build smoke test

Why:

- current CI only checks Ubuntu + Python 3.12
- packaging and PyO3 regressions often hide in version/platform gaps

Definition of done:

- `.github/workflows/test.yml` runs a matrix
- failures are visible per environment

### 2. Oracle enforcement (shipped 2026-04-19, reshaped from original plan)

Original plan (early 2026-04): a separate `oracle.yml` workflow on a
Mon/Wed/Fri cron to run `pytest tests/oracle/` + validation scripts.

What shipped: a workflow matching that plan landed in April and ran on
cron with `|| true` (advisory, not gating). Revisiting it on 2026-04-19
surfaced two problems:

1. Cron cadence + advisory mode meant regressions could sit on main
   for up to 2 days, silently. Removing `|| true` alone wasn't enough —
   the cadence was also wrong.
2. More importantly, running "oracle tests" as binary-pass-fail CI
   gates is a category error. Oracle tests are *measurements*, not
   assertions — a "failure" can be ferritin regressing, upstream
   drifting (BALL Julia's numbers shifted between April and 2026-04-18
   — caught by regeneration, not by reverting ferritin), a newly
   discovered convention gap, or tolerance miscalibration. None of
   those are things a green/red check can adjudicate.

What replaced it:

- **Library-parity tests** under `tests/oracle/` (ferritin-DSSP vs
  pydssp, ferritin-I/O vs biopython/gemmi, ferritin-SASA vs freesasa,
  ferritin-AMBER96 vs BALL-on-crambin) are really unit tests using
  another library as the ground-truth source. They now run on every
  PR via `.github/workflows/test.yml` (oracle.yml deleted).
- **Heavy oracle measurements** (OpenMM 1000-PDB AMBER96 parity,
  USAlign 4656-pair TM-align, 50K battle test, GPU kernels) live in
  `validation/`, run off-CI per release, and produce JSON/report
  artifacts. These are the scientific evidence — not gate-shaped.
- **When an oracle drifts**: the regeneration workflow is now a named
  procedure in `docs/ORACLE_SETUP.md §Regenerating reference values`,
  modelled on the 2026-04-19 BALL Julia regen.

The conceptual split (unit-test vs oracle-measurement) is captured in
the repo's oracle-testing philosophy (`devdocs/ORACLE.md`) and in the
reproducibility recipe (`docs/ORACLE_SETUP.md`).

### 3. Add Coverage Reporting

Track:

- Python line coverage
- Rust test counts at crate level

Initial requirement:

- report coverage first
- enforce thresholds later once baseline stabilizes

Why:

- test count alone hides untested surfaces
- current test distribution is strong but uneven

### 4. Add CLI Integration Tests

Add explicit tests for:

- `ferritin-ingest` single-file output
- `ferritin-ingest --per-structure`
- `tmalign` tabular output
- `usalign` tabular output
- failure modes on bad input

Why:

- ferritin’s data-wrangling story depends on the ingest and CLI layer
- these are product surfaces and should not rely on manual checks

---

## P1: Make Structure Preparation a Reliable Product Surface

### 5. Expose Public Prep APIs

Promote internal connector capabilities into the top-level Python package:

- `place_sidechain_hydrogens`
- `place_all_hydrogens`
- `place_general_hydrogens`
- `reconstruct_fragments`

Then add:

- `prepare_structure()`
- `batch_prepare_structures()`
- `load_and_prepare()`

Proposed pipeline:

1. load
2. reconstruct fragments
3. place hydrogens
4. minimize hydrogens
5. return structure + prep report

Prep report should include:

- atoms added
- residues skipped
- warnings
- whether ligands / waters / non-standard residues were touched
- energy before/after

Why:

- this is the missing interface for real use
- it aligns directly with the stated ferritin vision

### 6. Add End-to-End Prep Tests

Create integration tests for:

- complete standard protein prep
- missing-heavy-atom reconstruction
- idempotency of repeated prep
- prep on structures with ligands and waters
- prep on structures with non-standard residues
- tolerant behavior when partial prep is possible

Required assertions:

- no duplicate hydrogens
- atom counts change only when expected
- no NaN coordinates
- minimization reduces or stabilizes clash energy
- repeated prep does not keep mutating the structure

### 7. Add Prep Oracles

Compare against:

- BALL / BiochemicalAlgorithms.jl for hydrogen placement and reconstruction
- PDBFixer for full preparation behavior
- Gemmi for atom-count sanity where relevant

Focus first on:

- 20 standard amino acids
- disulfides
- termini
- waters
- common non-standard residues like MSE

---

## P2: Build a Real Regression Corpus

### 8. Create `tests/corpus/`

Add small, curated structures for:

- alt conformers
- insertion codes
- multi-model NMR
- missing sidechain atoms
- missing backbone atoms
- chain breaks
- ligands
- waters
- mixed protein/RNA if supported
- large-ish structures for performance smoke

Each fixture should have:

- a one-line reason it exists
- a linked bug or failure mode if applicable

### 9. Add Bug-Repro Tests

Policy:

- every parsing bug
- every alignment bug
- every prep bug
- every Arrow schema bug

gets a regression test before closing

### 10. Add Property / Fuzz Style Tests

Best candidates:

- selection parser
- Arrow schema conversion
- contact-map / distance-matrix invariants
- Kabsch / transform invariants
- load/save roundtrip on small generated structures

Good invariants:

- distance matrix symmetric, diagonal zero
- contact map monotonic with cutoff
- Kabsch RMSD never negative
- batch results match serial results
- roundtrip preserves counts and coordinates within tolerance

---

## P3: Strengthen Scientific Validation

### 11. SASA Validation

Keep current SASA path strong by formalizing:

- Bondi vs ProtOr comparisons
- FreeSASA library and CLI agreement
- tolerance by structure size class
- random corpus sampling with pinned seed

### 12. DSSP Validation

Current gap:

- DSSP has internal tests and sanity checks, but no strong external oracle in CI

Add:

- `mkdssp` or GROMACS-based oracle workflow
- residue-level agreement reports
- discrepancy categorization for chain ends, bends, and H-placement effects

### 13. H-Bond Validation

Split validation into:

- backbone H-bonds with energy criterion
- geometric H-bonds for broader atom classes

Compare against:

- DSSP where applicable
- GROMACS or mdtraj-style references if feasible

### 14. Forcefield / Minimization / MD Validation

Current status (2026-04-09):

- ForceField trait with AMBER96 and CHARMM19+EEF1 implementations
- Three minimizers: steepest descent, conjugate gradient, L-BFGS
- Full MD pipeline: Velocity Verlet, NVE/NVT, SHAKE/RATTLE, EEF1 solvation
- Torsion gradient, improper torsions, switching functions, neighbor list all implemented
- API behavior tested (28 Python tests + 44 Rust tests pass)
- Energy oracle comparison against BALL Julia available but not yet systematic

Add:

- energy component oracle comparison on small systems (BALL Julia on 50+ structures)
- L-BFGS vs CG vs SD convergence comparison on crambin
- long-run NVE drift tests (with SHAKE, with/without EEF1)
- CHARMM19 energy vs BALL Julia on shared test structures
- EEF1 solvation energy oracle vs BALL reference values
- minimization reproducibility tests
- known-clash relief test set

Be explicit:

- ferritin MD is suitable for structure relaxation and short equilibration
- for production MD, use GROMACS/OpenMM with ferritin for prep and analysis

---

## P4: Data-Wrangling Reliability

### 15. Make `ferritin-ingest` a Core Reliability Target

Add tests for:

- schema stability
- collision-safe structure IDs
- mixed PDB/mmCIF corpora
- partial failure handling
- chunked streaming writes
- roundtrip from Parquet back to structure where applicable

### 16. Add Per-Residue Feature Tables

Produce stable feature export for:

- residue SASA / RSA
- DSSP
- phi / psi / omega
- H-bond counts
- chain metadata
- optional interface labels

Why:

- this is the highest-value bridge from structures to downstream ML
- it reduces glue code and removes one of the main pain points for users

### 17. Version Schemas Explicitly

For Arrow/Parquet outputs:

- add schema version metadata
- test backward compatibility expectations
- define what constitutes a breaking schema change

---

## P5: Performance Reliability

### 18. Add Benchmarks to Protect the Fast Path

Track:

- load throughput
- SASA throughput
- alignment throughput
- ingest throughput
- memory on large corpora

Recommended approach:

- small stable benchmark set in repo
- optional larger benchmark on schedule

### 19. Prevent Known Scaling Traps

The roadmap already notes nonbonded `O(N^2)` forcefield behavior and missing neighbor lists.

Protect against regressions by:

- adding warning tests / documentation checks for large systems
- benchmarking pre- and post-neighbor-list implementation
- ensuring batch paths do not accidentally fall back to Python loops

---

## Concrete Deliverables

## Phase 1: Two Weeks

- add `RELIABILITY_ROADMAP.md`
- expand CI matrix
- add scheduled oracle workflow
- add CLI integration tests
- add coverage reporting

## Phase 2: Two to Four Weeks

- expose public prep APIs
- add end-to-end prep tests
- create regression corpus
- add schema/version tests for Arrow/Parquet

## Phase 3: One to Two Months

- add DSSP external oracle
- add BALL/PDBFixer prep comparisons on a broader corpus
- add per-residue feature export
- add performance regression reporting

---

## Suggested New Test Layout

```text
tests/
  test_alignment.py
  test_analysis.py
  test_arrow.py
  test_forcefield.py
  test_geometry.py
  test_hbond.py
  test_mmalign.py
  test_pdb_io.py
  test_sasa.py
  test_select.py
  test_prepare.py
  test_ingest_cli.py
  test_tmalign_cli.py
  test_usalign_cli.py
  oracle/
    test_io_oracle.py
    test_tmscore_oracle.py
    test_sasa_oracle.py
    test_dssp_oracle.py
    test_prepare_oracle.py
  corpus/
    README.md
    altloc/
    insertion_codes/
    missing_atoms/
    multimodel/
    ligands/
```

---

## Success Metrics

Ferritin should be considered reliability-ready when:

- every public API surface is covered by deterministic tests
- oracle validation runs automatically on a schedule
- CI covers supported Python versions and at least two OS families
- structure preparation has a public end-to-end workflow
- regression fixtures exist for major archive edge cases
- Arrow/Parquet schemas are versioned and tested
- benchmark regressions are visible, not anecdotal

---

## Recommended Immediate Next Actions

1. Update CI to add a Python/version matrix and a scheduled oracle workflow.
2. Add CLI integration tests for `ferritin-ingest`, `tmalign`, and `usalign`.
3. Expose structure-prep functions already present in the connector.
4. Implement `ferritin.prepare()` with a prep report object.
5. Build a fixed regression corpus from real problematic structures.

This order keeps ferritin focused on the actual user pain point:

**reliable, fast structure wrangling before downstream modeling or AI experimentation.**
