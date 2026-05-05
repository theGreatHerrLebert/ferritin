# Changelog

All notable changes to proteon are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The CHANGELOG is the **human narrative**; the machine-readable claims
manifest at `evident/evident.yaml` and per-release bundle at
`evident/reports/<tag>/manifest.json` are the **audit trail**. Each
release tag has a paired EVIDENT bundle pinned by sha256.

## [Unreleased]

## [0.1.4] — 2026-05-04

Second EVIDENT release. Trust pyramid completed: the dropped
fold-preservation claims (proteon CHARMM19+EEF1 vs OpenMM CHARMM36+OBC2,
proteon AMBER96 vs OpenMM AMBER96) are back, fully wired with
per-claim renderers. The 50K corpus oracle re-runs cleanly under
`pebble.ProcessPool` with **n_ok = 5,309 / 44,210** vs the v0.1.3
bundle's 807 — a 6.5× improvement on the same population. Manifest
trimmed from 20 → 18 force-of-evidence claims, then re-introduced the
two fold-preservation claims, landing at 20 release-tier claims with
real evidence.

### Added

- **Fold-preservation claims** (#50) — both proteon force fields are
  gated against OpenMM minimization on the same 1000-PDB random
  sample at TM-score level:
  - `proteon-charmm19-fold-preservation-vs-openmm-release-1k-pdbs`:
    proteon CHARMM19+EEF1 median TM=0.9944, OpenMM CHARMM36+OBC2
    median TM=0.9991, median tm_diff +0.0040, n_ok=886/1000.
  - `proteon-amber96-fold-preservation-vs-openmm-release-1k-pdbs`:
    proteon AMBER96 median TM=0.9959, OpenMM AMBER96 median
    TM=0.9992, median tm_diff +0.0028, n_ok=864/1000. Cleanest
    cross-implementation oracle proteon ships — both arms claim
    AMBER96, so the diff is purely implementation, not parameter
    set.
  - `validation/fold_preservation/join_fold_preservation.py` joins
    proteon-side and OpenMM-side per-PDB JSONLs into the canonical
    artifact `{pdb, n_ca, proteon: {...}, openmm: {...}, tm_diff,
    rmsd_diff_A}`.
  - `validation/report/render_fold_preservation.py` — per-claim HTML
    with TM distributions per side, tm_diff histogram against the
    claim's tolerance band, scatter, and top-20 outliers.
- **Three-tier capture schema** (#46) — `evident/CAPTURE_SCHEMA.md`
  defines `minimal` / `extended` / `full` capture levels for oracle
  runners, the per-PDB filesystem layout, the `hardware.json`
  schema, and the renderer contract. `PROTEON_CAPTURE_LEVEL` env
  controls the level.
- `validation/charmm19_eef1_ball_oracle.py` migrates from
  `concurrent.futures.ProcessPoolExecutor` to `pebble.ProcessPool`
  for true per-task subprocess isolation. Closes the
  `BrokenProcessPool` cascade that limited the v0.1.3 50K corpus
  run's `n_ok` to 807 of 44,210 attempted. (#44)
- New env knob `TASK_TIMEOUT_S` (default 60 s) on the corpus oracle
  runner — terminates hung BALL setup tasks individually instead of
  stalling the whole run.
- `pebble` added to both EVIDENT image variants' pip oracle install
  layer.

### Changed

- **Population definition for the 50K oracle is now "well-resolved
  wwPDB"** (#47). The runner skips PDBs with missing heavy atoms
  instead of running `PDBFixer.addMissingAtoms()`, which hangs
  deterministically on certain inputs. The defensible population
  shifts from "everything PDBFixer can repair" → "structures
  resolved well enough to score directly". 18,912 / 44,210 records
  on the v0.1.4 run skip for this reason — they are explicitly
  out-of-population, not failures.
- mkdssp build chain in both `evident/Dockerfile` and
  `evident/Dockerfile.cuda` is now version-pinned: `libmcfp v1.4.2`,
  `libcifpp v10.0.3`, `dssp v4.6.1`. Earlier `--depth 1 main` clones
  caused three breaking-upstream incidents during v1 prep. ARGs are
  exposed (`LIBMCFP_VERSION`, `LIBCIFPP_VERSION`, `DSSP_VERSION`)
  for downstream override. Closes #13.
- `evident/Dockerfile.cuda` switched from `jammy + 2 PPAs` to
  `noble + 0 PPAs` (#45) — eliminates the Launchpad/PPA flake class
  that broke 4 of the v0.1.3 cuda image builds.
- `libglib2.0-0` added to both runtime stages (#45). Apptainer smoke
  on monster3 caught `import ball` failing on
  `libglib-2.0.so.0: cannot open shared object file` in the cuda
  image; same fix applied to slim defensively.

### Removed

- 2 release-tier claims trimmed from `evident/evident.yaml` (#49)
  that asserted evidence we don't currently have:
  - `forcefield_amber_openmm` — removed; AMBER96 vs OpenMM AMBER96
    is now carried by the fold-preservation pair (#50).
  - `msa.yaml` — removed; MSA feature parity is held in unit-test
    scope, not framework-tier.

### Fixed

- 50K corpus run pass rate climbs from 807 / 44,210 (v0.1.3) to
  5,309 / 44,210 (v0.1.4) on the same input population. The lift
  is the sum of: pebble per-task isolation killing the
  `BrokenProcessPool` cascade, the `addMissingAtoms` hang fix,
  and the 60 s task timeout terminating the long tail of stuck
  BALL setups.

## [0.1.3] — 2026-05-03

First EVIDENT release: end-to-end claim / replay / report
architecture with sha256-pinned audit trail.

### Added

- **EVIDENT framework end-to-end** (Phases 1–4 across PRs
  #28 / #30 / #11 / #31):
  - `evident/scripts/lock_release_replays.py` walks every claim
    YAML under `evident/claims/`, sha256-pins each artifact, and
    emits an immutable bundle at `evident/reports/<tag>/`.
    Tier-aware: `pytest console output` style claims correctly
    surface as `ci-only` in the manifest, not `missing`.
  - Lock-time environment snapshot embedded in every manifest:
    Python version, platform, full `pip freeze`, sha256 of
    `cargo metadata`. Closes #36.
  - `evident/scripts/build_index.py` — multi-release aggregator
    served at `<pages-url>/evident/reports/`.
  - `evident/scripts/cut_release.sh` — one-command flow: scp
    artifacts from monster3 → lock bundle → print git tag steps.
  - Per-claim HTML renderers for the corpus oracle, SASA-vs-Biopython
    1k, and 50K battle test. Reports embed plots as base64,
    self-contained, no CDN dependency.
- **Reproducibility images** on GHCR
  (`ghcr.io/thegreatherrlebert/proteon-evident:v0.1.3`):
  - `:slim` — Debian trixie + Python 3.13. Bundles proteon,
    `ball-py`, `openmm`, `pdbfixer`, `biopython`, `gemmi`,
    `freesasa`, `pydssp`, `gromacs`, `mmseqs2`, `mkdssp`,
    `USalign`, `reduce`.
  - `:cuda` — nvidia/cuda 12.8.2 + Ubuntu noble + Python 3.12.
    Same payload, GPU-enabled `proteon-connector` + `cudarc`.
  - Entrypoint dispatches `replay <claim-id>`,
    `render <release-tag>`, plus pass-through to the vendored
    `evident` CLI.
- **CHARMM19+EEF1 oracle** (PRs #16–#26):
  - All seven force-field components (bond_stretch, angle_bend,
    proper_torsion, improper_torsion, vdw, electrostatic,
    EEF1 solvation) gated against BALL on crambin within
    per-component bands.
  - 50k-ready corpus runner with chunked-pool segfault isolation
    (later replaced by pebble in v0.1.4) and a protein-only
    filter that drops nucleic-acid PDBs the force field can't
    parametrise.
  - Distance-dependent dielectric (ε ∝ r) — the canonical
    CHARMM19 convention — replaces the prior constant-ε
    implementation.
  - PHE/TYR aromatic-ring para-diagonal LJ exclusion port from
    BALL's `charmmNonBonded.C:547-565`.
- **CHARMM19+EEF1 corpus oracle**:
  - 1k-PDB BALL oracle on the curated `validation/pdbs/`. Released
    alongside the unit oracle as `proteon-charmm19-vs-ball-corpus-1k-pdbs`.
  - 50K-PDB BALL oracle on a random wwPDB sample, run on monster3.
    Headline `n_ok=807 of 44 210 attempted` reflects the
    `BrokenProcessPool` cascade; per-component bands on the 807
    successful records match the unit oracle (median 0.05%
    bond_stretch, 0.001% torsion, 0.02% electrostatic, 1.62%
    EEF1 solvation, 0.41% improper). Cascade noise framed
    honestly in the claim's `failure_modes`; pebble migration
    is the v0.1.4 follow-up.
- **EVIDENT manifest** (`evident/evident.yaml`) wires 20 claims
  across all proteon subsystems (force fields, alignment, SASA,
  DSSP, hydrogens, GB OBC, supervision, MSA, I/O, pipeline batch).
  Each claim names its oracle, tolerance, replay command, and
  artifact path.

### Changed

- 11 of 13 CI-tier oracle claims now pass against external
  references (BALL, Biopython, USAlign, pydssp, reduce). The two
  pending — `proteon-amber96-vs-openmm-release-1k-pdbs` and
  `proteon-msa-vs-mmseqs2-research` — have monster3-side artifacts
  not yet mirrored.

### Removed

- 6 low-signal claims trimmed in PR #43 (`acceleration_paths`,
  `gpu_cpu_parity`, `forcefield_charmm19_internal`,
  `sasa_freesasa`, both fold-preservation entries). The two
  fold-preservation claims drop until #37 (PDBFixer pre-pass on
  the GROMACS runner) lands real artifacts. Manifest goes from
  26 → 20 claims; every "missing" row is now a real coverage gap,
  not a documentation artefact.

### Known gaps tracked for follow-up

- **#37** — fold-preservation GROMACS artifact regen needs
  PDBFixer pre-pass + monster3 GROMACS install.
- **#42** — v0.2.0 vision: every release-tier claim at 50K scale
  across all oracles.
- **#44** — pebble runner migration (resolved in [Unreleased]).

## [0.1.2] — 2026-04-24

First clean-room MIT release after the 2026-04-24 third-party
lineage audit.

### Changed

- Phase 3 hydrogen placer in `proteon-connector/src/add_hydrogens.rs`
  rewritten clean-room from standard crystallographic geometry; the
  prior dispatcher was structurally derived from BALL's
  `AddHydrogenProcessor` (LGPL-2.1).
- `reconstruct.rs` reattributed to its actual upstream: MIT
  `BiochemicalAlgorithms.jl`, not LGPL BALL C++.
- `forcefield/md.rs` and `forcefield/gb_obc.rs` cite primary
  literature / OpenMM's MIT Reference Platform.

### Added

- `THIRD_PARTY_NOTICES.md` consolidates verbatim upstream license
  texts for TM-align/US-align, MMseqs2, BiochemicalAlgorithms.jl,
  OpenMM, and pdbtbx.

### Performance

- SASA GPU auto-dispatch threshold raised from 500 → 10 000 atoms.
  Small-protein batch throughput recovered ~10× against the
  pre-2026-04-12 regression (monster3 5K: 18.6/s → 192.5/s; numerical
  output bit-identical).

### Benchmark

- `benchmark/run_benchmark.py` now tracks `n_negative_energy`
  alongside `n_converged`. After the 2026-04-11 CHARMM19 heavy-atom
  relaxation default, negative final energy is the documented
  correctness invariant; `converged` is preserved for back-compat.

### No public API changes.

## [0.1.1] — 2026-04-18

First clean PyPI release with populated package metadata.

### Added

- `pip install proteon` installs the Pythonic wrapper; the
  PyO3 `proteon-connector` is pulled in transitively.
- Wheels published for Linux x86_64 / aarch64, macOS x86_64 /
  aarch64, and Windows x86_64.

### Validation highlights at the time of release

- AMBER96 matches OpenMM to within 0.2% on every energy component
  at NoCutoff.
- CHARMM19+EEF1 pipeline: 50K random PDBs at 99.1% end-to-end
  success in 3.5h on RTX 5090.
- TM-align port: 0.003 median TM-score drift from the reference
  C++ USalign across 4 656 pairs.
- SASA: 0.17% median deviation vs Biopython on 1 000 structures.

### No runtime code changes since 0.1.0.

## [0.1.0]

Initial PyPI publication (one-shot to claim the project name; pages
rendered with empty descriptions, superseded by 0.1.1).

[Unreleased]: https://github.com/theGreatHerrLebert/proteon/compare/v0.1.3...HEAD
[0.1.3]: https://github.com/theGreatHerrLebert/proteon/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/theGreatHerrLebert/proteon/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/theGreatHerrLebert/proteon/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/theGreatHerrLebert/proteon/releases/tag/v0.1.0
