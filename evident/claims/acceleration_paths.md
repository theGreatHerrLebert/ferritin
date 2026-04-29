# Acceleration-path parity: SIMD, rayon, multi-thread supervision (gap)

Operational case writeup for the acceleration-paths-parity
*reference* claim in `claims/acceleration_paths.yaml`.

## Why this is a gap, not a measurement

`feedback_cross_path_parity.md` (memory) states the principle:

> Any accelerated code path (NBL, SIMD, GPU, rayon) needs a parity
> test against the slow path on every component, parametrized over
> force field.

Today's coverage:

- **NBL parity** ✓ (`tests/test_cross_path_parity.py`, claim
  `proteon-cross-path-parity-nbl-ci`).
- **GPU parity** — separate gap claim
  (`proteon-gpu-vs-cpu-energy-parity-reference`).
- **SIMD parity** — flagged in cross-path-parity failure_modes;
  no asserting test.
- **rayon vs serial parity** for energy compute — same.
- **multi-thread parity** for the supervision Rust extractor —
  flagged in `proteon-supervision-rust-vs-python-ci` failure_modes;
  current test runs single-thread only.

This claim is the placeholder for the three remaining axes: SIMD,
rayon, and multi-thread supervision. They share enough structure
to bundle: each is a parametrisation of an existing parity test
plus a code-path forcing flag.

## Required test design

Three sub-tests, all kind=reference until they land:

### SIMD parity
Extend `tests/test_cross_path_parity.py` (or add
`tests/test_simd_parity.py`) to force scalar-only computation
(via a build flag or env var that disables SIMD intrinsics) and
compare against the SIMD-enabled path on the same 4-structure × 2-FF
matrix. Tolerance: same `max(1e-6 kcal/mol, 1e-9 * |E|)` as NBL
parity.

### rayon vs serial parity
Same fixture matrix, run with `n_threads=1` (forces serial) and
`n_threads=-1` (forces all cores). Each component must agree to
the same tolerance. Note proteon-internal trap: `n_threads=0` runs
*serial*, not all-cores (per memory `feedback_proteon_n_threads_trap`).
The test must use `1` and `-1`, not `0`.

### Multi-thread supervision parity
Extend `tests/test_supervision_rust_parity.py` to parametrise over
`n_threads` ∈ {1, -1}. The Rust supervision extractor is
data-race-free but reduction-order may be non-deterministic under
rayon; tolerance same as the existing supervision claim
(`assert_array_equal` on coords/masks, `atol=1e-5` on angles/frames).

## Trust strategy

Same Reference Shadowing pattern as the existing cross-path-parity
and supervision-rust-vs-python claims. The slow / serial / single-
thread path is the reference; the accelerated path is the
optimised one.

## What this claim does NOT say

- It does not assert that any of these three accelerated paths is
  currently broken. The 50K battle test, kernel-internal Rust
  tests, and production usage all give *some* evidence of
  correctness — just not the structured public-API parity
  evidence.
- It does not block release. The CPU/serial/single-thread paths
  are exercised by their own asserting claims; users worried
  about the acceleration paths can disable them via build flags
  or thread budgets.

## How to close

Three small PRs, in priority order:

1. **rayon parity** — easiest, cheapest, highest leverage (rayon
   is on by default). Promote to its own
   `proteon-rayon-vs-serial-parity-ci` measurement claim.
2. **multi-thread supervision parity** — second, since the
   supervision pipeline runs at production scale during corpus
   builds. Promote to its own measurement claim.
3. **SIMD parity** — last, because SIMD is gated by build flags
   that most users don't toggle. Still worth closing.

Each closure converts this bundled reference claim into one
measurement claim per axis, and the failure_modes notes in
existing claims (cross-path-parity, supervision-rust-vs-python)
get downgraded to "covered in sibling claim".
