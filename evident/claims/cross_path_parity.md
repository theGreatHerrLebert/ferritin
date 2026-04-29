# Cross-path parity: NBL vs exact

Operational case writeup for the cross-path-parity CI claim in
`claims/cross_path_parity.yaml`.

## Problem

Proteon's `compute_energy` has two implementations of the same math:
- An **exact** O(N²) all-pairs path (`compute_energy_impl`), used for
  small structures and forced via `nbl_threshold=10**9`.
- A **neighbour-list (NBL)** path (`compute_energy_nbl`), used for
  larger structures (default heuristic threshold = 2000 atoms) and
  forced via `nbl_threshold=0`.

A regression in either path that shifts a per-component energy
without crossing any invariant threshold (still finite, still
non-zero, still has the right sign) would ship invisibly on a
single-path test. **Cross-path parity** — every component, every
force field, every reference structure — is the highest-leverage
guard against that class of bug.

The 2026-04-11 CHARMM+EEF1 history is the worst-case demonstration:
`compute_energy_and_forces_nbl()` was silently skipping the EEF1
solvation term entirely while `compute_energy_impl` returned the
real value. The divergence was O(10³) kcal/mol and would have failed
this test on day one, but the bug shipped because no parity test
existed at the time. This claim is the regression guard that closes
the door.

## Trust Strategy

Validation against the same code's slow reference path. This is *not*
an external-oracle claim — the "oracle" is proteon's own
exact O(N²) implementation, treated as the reference for the
accelerated NBL path. That makes it a Reference Shadowing pattern in
EVIDENT terms: maintain a simple reference implementation alongside
the optimised version and assert agreement.

- **Reference path**: `nbl_threshold=10**9` forces
  `compute_energy_impl` (exact O(N²) all-pairs).
- **Accelerated path**: `nbl_threshold=0` forces
  `compute_energy_nbl` (neighbour-list).
- **Default path**: no `nbl_threshold` argument; selects exact below
  2000 atoms, NBL above.
- **Corpus**: 4 reference structures spanning the threshold —
  1crn (327 atoms), 1bpi (454), 1ubq (602), 1ake (3317). The
  one-large structure (1ake) exercises the NBL path under the
  default heuristic; the three small ones exercise the exact path.
- **Force fields**: amber96, charmm19_eef1. Parametrised over both
  → 4 × 2 = 8 measurement cells per parity assertion.

## Evidence

`tests/test_cross_path_parity.py`:

- `test_nbl_matches_exact_all_components` — for each (structure, FF),
  computes `e_exact` and `e_nbl` and asserts on every component
  (`bond_stretch`, `angle_bend`, `torsion`, `improper_torsion`, `vdw`,
  `electrostatic`, `solvation`, `total`) that
  `|exact - nbl| <= max(TOL_KJ_MOL, 1e-9 * |exact|)` where
  `TOL_KJ_MOL = 1e-6 * 4.184` (1e-6 kcal/mol in kJ/mol). The
  relative term scales the floor for large electrostatic components
  on clashy raw PDBs (O(10⁵) kJ/mol).
- `test_default_path_matches_one_of_the_forced_paths` — the default
  (no override) must match `force_nbl` exactly for >2000 atoms and
  `forbid_nbl` exactly for ≤2000 atoms, to the same tolerance. This
  pins the heuristic itself as a contract.

The Rust side has a parallel parity test in
`proteon-connector/src/forcefield/energy.rs` gradient_tests; the
Python test exercises the full PyO3 binding path including kJ/mol
unit conversion and dict assembly, so both sides close the loop from
public API down to kernel.

## Assumptions

- `1e-6 kcal/mol` is the right strictness floor: NBL vs exact drift
  is either zero (correct) or O(10⁰+) kcal/mol (a bug); there is no
  physical process producing drift in between. Tighter would be
  flaky on FP reassociation; looser would mask real bugs.
- The relative `1e-9 * |exact|` term scales for large-magnitude
  components on raw PDBs without minimisation; minimised structures
  rarely need it.
- The 4-structure corpus spans both regimes of the default heuristic
  (≤ and > 2000 atoms). A heuristic change (e.g. lowering to 1000)
  would break `test_default_path_matches_one_of_the_forced_paths`
  on at least one structure, surfacing it as a manifest-visible
  regression.
- `compute_energy_impl` is the correct reference. It is the simpler
  implementation (no neighbour list, no chunking, no SIMD); a bug
  in it would invalidate the parity test direction. The invariant
  test suite (sign checks, finiteness, total = sum components)
  guards `compute_energy_impl` independently.

## Failure Modes

- **Symmetric bugs invisible.** If both paths share the same wrong
  formula (the 2026-04-11 EEF1 1-2/1-3 exclusion bug), the parity
  test passes — both paths agree on the wrong answer. The invariant
  suite's sign and magnitude checks plus the BALL/OpenMM oracle
  claims close that gap.
- **No GPU path here.** This claim covers only NBL-vs-exact CPU
  parity. GPU parity (CPU to 1e-11 per memory) is a separate
  cross-path claim that should be added when the GPU path's tests
  are wired into a Python parity check at the public API.
- **No SIMD/rayon path here.** Same as above — those acceleration
  paths need their own parity claims, parametrised over the same
  4-structure × 2-FF matrix.
- **Heuristic threshold is a contract.** `NBL_AUTO_THRESHOLD = 2000`
  is asserted via `test_default_path_matches_one_of_the_forced_paths`
  but only at the boundary of the corpus (≤602 vs 3317). A threshold
  change to e.g. 1500 would not surface unless a new fixture lands
  between 1500 and 2000 atoms.
- **Solvation ambiguity for AMBER96.** AMBER96 has zero solvation; the
  parity check still runs (both paths return 0.0), so a regression
  that introduced a non-zero AMBER96 solvation would surface as
  divergence between paths and as a sign/magnitude invariant
  failure.

## Lessons

- Reference Shadowing pays for itself the moment an acceleration
  path diverges. The 2026-04-11 EEF1 incident is the canonical
  story: a silent O(10³) kcal/mol divergence that would have been
  caught immediately by this exact test.
- Parametrise over the full FF × structure × path matrix once, not
  per-test. The CHARMM and AMBER paths share enough infrastructure
  that a single parity test parametrised over both is more
  defensible than two separate ones.
