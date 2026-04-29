# Supervision: Rust batch extractor vs Python reference

Operational case writeup for the supervision Rust-vs-Python parity
CI claim in `claims/supervision_parity.yaml`.

## Problem

Proteon's Layer 5 geometric-DL data export has two implementations
of the same supervision schema:

- **Python reference**: `proteon.batch_build_structure_supervision_examples`
  (lives in `supervision.py` / `supervision_geometry.py`). The
  authoritative reference; what the schema *means* operationally.
- **Rust hot path**:
  `proteon_connector.py_supervision.batch_extract_structure_supervision`.
  An accelerated batch extractor used by default when the Rust
  backend imports cleanly; falls back to the Python reference
  otherwise.

A divergence between the two paths silently poisons downstream
training: half a corpus might be built under Rust (fast, but
regressed) and half under Python (correct), with the model learning
the Rust artifacts as features. The trust question is therefore
**bit-for-bit parity** on every output tensor across both paths,
on a small fixture set that exercises every supervision computation
class:

- Coordinate extraction (all_atom_positions, atom14_gt_positions,
  pseudo_beta).
- Mask computation (atom14_gt_exists, pseudo_beta_mask, plus per-
  angle masks).
- Torsion math (phi/psi/omega/chi — all four dihedral kinds).
- Rigidgroup frame construction (8 groups per residue).
- Batch padding semantics (Rust pads to longest in batch; Python
  returns each example at its real length).

## Trust Strategy

Validation by Reference Shadowing on the same code's slow reference
path. The Python implementation is the reference; the Rust
implementation is the optimised one. Same EVIDENT shape as the
cross-path-parity claim, just on tensor outputs instead of energies.

- **Reference path**: Python in `proteon.supervision*`.
- **Accelerated path**: Rust in `proteon_connector.py_supervision`.
- **Corpus**: 1crn (46 res, real sidechains — exercises chi-angle
  + rigidgroup-frame paths) and 1ubq (76 res, broader residue
  composition for chi-angle atom-lookup coverage).
- **Batch test**: 1crn + 1ubq mixed-length through the batch path,
  asserts Rust's pad-then-slice matches Python's unpadded per-
  example output.

The test is `pytest.mark.skipif(not rust_supervision_available())` —
on a CPU-only Python environment without the connector, the test
skips cleanly and the Python path is the only path. The Rust parity
matters specifically when both are present (any production
environment).

## Evidence

`tests/test_supervision_rust_parity.py::TestRustPythonParity` runs
three parametrised tests:

- `test_crambin_rust_matches_python` — 1crn through both paths.
- `test_ubiquitin_rust_matches_python` — 1ubq through both paths.
- `test_mixed_length_batch_rust_matches_python` — both PDBs through
  the batch path, slicing each example to its real length before
  comparison.

For each (structure, field group) cell:

- **Coordinate fields** (`all_atom_positions`, `atom14_gt_positions`,
  `pseudo_beta`): `np.testing.assert_array_equal` — bit-for-bit.
- **Mask fields** (10 total — `all_atom_mask`, `atom14_gt_exists`,
  `atom14_atom_exists`, `atom37_atom_exists`, `pseudo_beta_mask`,
  `rigidgroups_gt_exists`, `phi_mask`, `psi_mask`, `omega_mask`,
  `chi_mask`): byte-exact via `assert_array_equal`.
- **Angle fields** (`phi`, `psi`, `omega`, `chi_angles`):
  `assert_allclose(atol=1e-5, rtol=1e-5)`. Allowing a tiny float
  tolerance because dihedral math involves sqrts and atan2; in
  practice (2026-04-14 crambin + ubiquitin) diff == 0.
- **Frame fields** (`rigidgroups_gt_frames`):
  `assert_allclose(atol=1e-5, rtol=1e-5)`. Same rationale.

## Assumptions

- `assert_array_equal` is the right strictness for coords and masks:
  the supervision pipeline is deterministic on the same input atoms,
  so any difference is a real bug.
- The 1e-5 tolerance on angles and frames absorbs Rust-vs-Python
  sum reordering in the dihedral math without permitting drift.
  Empirical diff is 0 on crambin and ubiquitin; the floor exists
  in case a future residue triggers a path where one side reorders.
- Crambin + ubiquitin together cover the chi-angle atom-lookup
  paths for the 20 standard amino acids (THR/CYS/PRO/VAL/LEU in
  crambin; ubiquitin adds the rest of the standard composition).
  A residue not present in either fixture (selenomethionine,
  pyrrolysine, modified residues) would not be exercised.
- The 8 rigidgroups per residue are exercised on both fixtures;
  no fixture-specific rigidgroup gap.

## Failure Modes

- **Skip-on-missing-connector blind spot.** The test
  `pytest.mark.skipif`s when `rust_supervision_available()` is
  False. CI must positively assert that the Rust connector is
  importable, or a CPU-only environment can pass the test suite
  while exercising none of this claim. A regression in the build
  that hides the Rust path would silently green-light the Python
  fallback.
- **Two structures only.** No exotic-residue coverage (Sec, Pyl,
  modified amino acids). A regression in chi-angle lookup specific
  to those would not surface.
- **Bit-exact masks could become brittle.** Masks are integer-
  typed; the strict equality is appropriate. But a future change
  to the supervision schema (e.g. adding an "asymmetric unit"
  mask) would need both sides updated in lockstep, with this
  parity test as the gate.
- **Batch padding semantics.** `test_mixed_length_batch` verifies
  that Rust's pad-then-slice matches Python's unpadded output
  inside the valid-length region. A bug in how zeros bleed into
  the valid region (off-by-one slice) would surface; a bug
  outside the valid region would not (we only compare to
  `:py_example.length`).
- **Single-batch, single-process.** Tests run on a single
  process. Rust tensors that look correct on a single-threaded
  call but drift under rayon parallelism (data-race-free but
  reduction-order-dependent) would not surface here. A separate
  test running the batch extractor with `n_threads > 1` would
  close that gap.

## Lessons

- Reference Shadowing on tensor outputs has the same load-bearing
  role as on scalar energies. The cross-path-parity claim guards
  the energy compute; this claim guards the supervision compute.
  Together they pin the two acceleration paths most likely to
  drift silently.
- Bit-exact for coords and masks, near-exact for angles and
  frames is the right tolerance shape. Loosening coords/masks to
  match the angle tolerance would let drift hide; tightening
  angles to bit-exact would be flaky on FP reassociation across
  Rust and Python sum loops.
