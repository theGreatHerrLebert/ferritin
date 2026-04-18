# AMBER_UPDATE — 2026-04-14

Session-end snapshot of the AMBER96 + OBC GB investigation, plus the
surrounding Layer 5 / GPU work that landed alongside it.

---

## The investigation: heavy-atom displacement gap

The validation HTML report (`validation/report/report.html`) called out a
clear divergence in heavy-atom displacement between ferritin and the
other pipelines:

- Ferritin CHARMM19+EEF1: median RMSD **0.585 Å**
- OpenMM CHARMM36+OBC2:    median RMSD **0.213 Å**
- Ferritin AMBER96:        median RMSD **0.497 Å**
- OpenMM AMBER96+OBC:      median RMSD **0.218 Å**
- GROMACS AMBER96:         median RMSD **0.170 Å**

The report's current explanation (polar-H united-atom absorbing
implicit hydrogens) is **CHARMM-specific** and doesn't transfer to
AMBER96, where both sides are all-atom.

### What I confirmed

Per-residue diagnostic on clean X-ray fixtures (1crn, 1ubq, 1bpi,
1aaj — `constrain_heavy=False`, 100 minimization steps to match the
1000-PDB benchmark config):

| PDB  | L   | AMBER96 heavy RMSD | CHARMM heavy RMSD |
|------|-----|--------------------|-------------------|
| 1crn | 46  | 0.036 Å            | 0.728 Å           |
| 1ubq | 76  | 0.036 Å            | 0.763 Å           |
| 1bpi | 58  | 0.034 Å            | 0.804 Å           |
| 1aaj | 105 | 0.034 Å            | 0.020 Å           |

On clean X-ray inputs ferritin AMBER96 barely moves anything —
0.03–0.04 Å is thermal noise. The 0.497 Å median in the 1000-PDB
benchmark comes from structures with strained starting geometry
(reconstructed atoms, NMR sets, low-res crystals).

### What's actually causing the gap

For the **AMBER96 comparison**, the dominant driver is **missing
implicit solvent**:

- OpenMM benchmark: AMBER96 + **OBC GB** (`amber96_obc.xml`)
- Ferritin: AMBER96 in **vacuum**

OBC screens Coulomb forces. With GB the force landscape is gentler,
so OpenMM's LBFGS converges with smaller atom movement. In vacuum,
strained inputs have larger residual gradients and atoms move further.

For the **CHARMM comparison**, the gap is compound: different FF
family (19 polar-H vs 36 all-atom) AND different solvent (EEF1 vs
OBC2). Two confounded differences.

---

## OBC GB Phase A: shipped today

Commit `473f7cf`. Module + parameters + trait hooks + a calibrated
failing oracle. The pair-integral math, Born-radii rescaling, and
analytical forces are deliberately `todo!()` pending a dedicated
math session.

### Why scaffold-first

GB analytical forces (chain rule through per-pair Born-radii
dependence) are where sign errors slip in. Landing the contract +
calibrated failing oracle first forces the math session to converge
against a verified spec rather than "looks right" heuristics.
Ferritin's value-prop is "oracle-validated"; shipping plausible-
but-unchecked GB physics would undo that.

### What's in

`ferritin-connector/src/forcefield/gb_obc.rs`:
- `ObcAtomParams { radius, scale }` — per-atom Bondi-style radius +
  HCT overlap scale factor.
- `ObcGbParams::obc1()` / `obc2()` — pinned constants:
  - OBC1 (α=0.8, β=0.0, γ=2.909125) matches OpenMM `amber96_obc.xml`
  - OBC2 (α=1.0, β=0.8, γ=4.85) matches `charmm36_obc2.xml`
  - ε_in=1.0, ε_out=78.5, offset=0.09 Å on both
- Function signatures mirroring the EEF1 entry pattern in `energy.rs`:
  `compute_born_radii`, `gb_obc_energy`, `gb_obc_energy_and_forces`,
  `gb_obc_energy_and_forces_nbl`. Bodies are stubs that compile
  cleanly so dispatch wiring (Phase C) can be tested end-to-end
  before the math lands.
- Doc-comment carries the full OBC1 derivation: HCT integral,
  ψ = I·ρ', R_eff = 1/(1/ρ' − tanh(αψ − βψ² + γψ³)/ρ),
  GB pair formula f_GB = √(r² + R_iR_j exp(−r²/(4R_iR_j))).

`ferritin-connector/src/forcefield/params.rs`:
- `ForceField` trait gains `get_obc_gb()` / `has_obc_gb()`,
  default-False so existing CHARMM/AMBER impls compile unchanged.

`validation/amber96_obc_oracle.py`:
- PDBFixer-add-H → temp PDB → both pipelines load the identical
  hydrogens (only variable is each tool's energy implementation).
- OpenMM: `amber96.xml + amber96_obc.xml`, NoCutoff, two contexts
  (with + without GBSAOBCForce) → vacuum vs total decomposition.
- Ferritin: `compute_energy(ff="amber96", nbl_threshold=1e9)` for
  exact O(N²) NoCutoff equivalent.
- Tolerance gate: `|Δtotal/E_om| < 1%` AND `|ΔGB/GB_om| < 5%`.

Per-atom defaults documented in `gb_obc.rs` (Phase B loads from
`amber96_obc.ini`):
- H: radius=1.2 Å, scale=0.85
- C: radius=1.7 Å, scale=0.72
- N: radius=1.55 Å, scale=0.79
- O: radius=1.5 Å, scale=0.85
- S: radius=1.8 Å, scale=0.96
- P: radius=1.85 Å, scale=0.86

### Calibrated baseline on crambin (this is the Phase D target)

```
OpenMM AMBER96+OBC1 NoCutoff:
  total:  7982.734 kJ/mol
  vacuum: 9255.723 kJ/mol
  GB:    -1272.989 kJ/mol  (1 GB force removed)

Ferritin AMBER96 (vacuum, GB stub):
  total:  9179.122 kJ/mol
  vacuum: 9179.122 kJ/mol
  GB:        0.000 kJ/mol

Comparison:
  vacuum Δ:   76.6 kJ/mol  (0.83% — within existing AMBER96 oracle)
  total Δ:  1196.4 kJ/mol  (14.99% — entirely the missing GB term)
  GB Δ:     1273.0 kJ/mol  (100% — exact stub-vs-real gap)
```

The vacuum side already tracks within 0.83%, so AMBER96's existing
parity is preserved. The full 15% gap collapses to the missing GB
term. When Phase B/C land and `compute_born_radii` + `gb_obc_energy`
return −1273 kJ/mol on this molecule, total Δ should drop to <1% and
the script exits with `PHASE D PASS`.

### Tests

- 3 Rust unit tests pin the OBC1/OBC2 constants and the τ formula —
  any future edit that drifts them fails immediately.
- All 80+1 existing forcefield tests still pass.

---

## What remains (Phase B → E)

Estimated 6–8h of focused work split across phases.

### Phase B: math implementation (~3–4h)

1. **HCT pair integral + OBC1 tanh rescaling** in `compute_born_radii`.
   Reference: OpenMM `platforms/reference/src/SimTKReference/ReferenceObc.cpp`.
   Verify by computing Born radii on a few small systems and comparing
   to OpenMM single-point.
2. **`gb_obc_energy`** (E_self + E_pair via f_GB). After this lands,
   `validation/amber96_obc_oracle.py` should match within ~1% if
   Born radii are right.
3. **Analytical forces in `gb_obc_energy_and_forces`**. Verify against
   finite-difference on every atom for crambin before declaring done —
   sign errors here are the dominant failure mode.

### Phase C: dispatch wiring (~30 min)

Mirror the EEF1 hook pattern in `energy.rs`. Three call sites: lines
233 (energy-only), 482 (energy+forces), 612 (NBL). Route through
`result.solvation` so existing downstream code sees it.

### Phase D: oracle + parity (~30 min)

Run `amber96_obc_oracle.py` to `PHASE D PASS`. Add cross-path NBL
parity test (per the existing `feedback_cross_path_parity` memory:
any accelerated path needs a parity test against the slow path on
every component).

### Phase E: benchmark rerun + report fix (~30 min)

Re-run `validation/tm_fold_preservation_amber.py` on monster3 with
GB-enabled ferritin AMBER96. Expect median RMSD to drop from
0.497 Å toward OpenMM's 0.218 Å, closing the gap that prompted this
investigation.

Update `validation/report/report.html` narrative — the current
"polar-H united-atom" explanation is CHARMM-only and needs a separate
AMBER paragraph (the actual driver was the missing OBC GB, not H
absorption).

---

## Surrounding work that landed today

Context for whoever picks up Phase B+:

### GPU search (warp-collab PSSM Smith-Waterman)

- Phase 4.5a (`825d72e`) singletile warp-collab kernel — 9.67× on 5090
  for short queries (≤256 residues).
- Phase 4.5b (`55212f8`) multitile variant — 124× on 5090 for queries
  300–800, 74× for queries 500–1500.
- Covers the entire typical-protein-query range on GPU.

### Layer 5 (geometric DL data pipeline)

- Baseline `ac3b347` (was sitting uncommitted with ~1400 LOC + 540 LOC tests)
- Tensor SHA-256 checksums (`33a6226`)
- 10-class failure taxonomy (`58f9f1e`)
- Real-PDB smoke (`2d1f7d5`)
- Rust-Python supervision parity guard (`97fc977`)
- Three review findings fixed (`7713503`)
- GPU MSA search wired into SequenceExample (`4536a4f`)
- Sequence-based template features (`8cc6d0a`)
- Training NPZ export with checksum (`25db0a5`)
- First real end-to-end search release driver (`b7a5fc2`)

**69/69 Python Layer 5 tests, 173/173 ferritin-search Rust+CUDA tests,
80+3 forcefield Rust tests.**

### Demonstrated end-to-end

`validation/first_real_search_release.py`: 4 PDBs (1crn, 1ubq, 1bpi,
1aaj) through the full stack in 12.6s. MSA cross-hits, templates with
self-hit exclusion, training NPZ checksummed, 0 validation issues,
0 ingestion failures.

---

## Files of interest for next session

- `ferritin-connector/src/forcefield/gb_obc.rs` — Phase A scaffold
  with the algorithm spec in the doc-comment
- `ferritin-connector/src/forcefield/energy.rs:733-918` — EEF1
  reference implementation to mirror
- `validation/amber96_obc_oracle.py` — calibrated forcing function
- `validation/amber96_oracle.py` — existing AMBER96 vacuum oracle
  pattern (load via PDBFixer, write temp PDB, both pipelines read it)
- OpenMM `platforms/reference/src/SimTKReference/ReferenceObc.cpp`
  (in your OpenMM source tree) — canonical math reference
- Onufriev, Bashford, Case (2004) — the OBC paper
