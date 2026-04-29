# GPU vs CPU energy parity (gap)

Operational case writeup for the GPU↔CPU energy parity *reference*
claim in `claims/gpu_cpu_parity.yaml`.

## Why this is a gap, not a measurement

Proteon ships a `cuda` feature that auto-dispatches CHARMM19+EEF1
energy + SASA to GPU when a usable device is detected
(`proteon-connector/src/forcefield/gpu.rs`). Memory cites
"GPU matches CPU to 1e-11 on crambin" (OBC GB Phase B,
2026-04-15) and the 50K production battle test ran on the GPU code
path with a 99.1% physics-consistency rate.

But: **no Python-level test asserts GPU↔CPU parity at the public
API**. The Rust side has parity tests in
`proteon-connector/src/forcefield/gpu.rs` test modules; those
exercise kernel-internal correctness on synthetic inputs, not the
end-to-end `compute_energy(..., ff="charmm19_eef1")` Python path.
Users running production GPU pipelines have no manifest-visible
parity assertion to point at.

This claim is a `kind: reference` placeholder marking the gap. It
becomes a `kind: measurement` claim the day a Python parity test
lands.

## Required test design

A new test, parametrised over the same fixture × FF matrix as
`tests/test_cross_path_parity.py`, that:

1. Detects whether `cuda` feature is built and a device is
   available (skip cleanly otherwise — same pattern as existing
   `pytest.importorskip` and `REDUCE_BIN` env gates).
2. Forces the GPU code path via `proteon.compute_energy(...,
   nbl_threshold=0)` on a `cuda`-built proteon, and records every
   energy component.
3. Forces the CPU exact path on the same input and same
   parameters.
4. Asserts `|gpu - cpu| <= max(TOL_KJ_MOL_GPU, 1e-9 * |cpu|)` per
   component, where `TOL_KJ_MOL_GPU` is calibrated against the
   memory-cited 1e-11 baseline plus a margin for FP reassociation
   on different hardware.

The fixture set should span the same atom-count regimes as
cross-path-parity (1crn 327, 1bpi 454, 1ubq 602, 1ake 3317) — and
should additionally cover at least one structure that exceeds GPU
memory for a single chain to exercise the chunking path.

## Trust strategy

Validation by Reference Shadowing — the CPU exact path is the
reference, the GPU path is the optimised one. Same shape as the
existing `proteon-cross-path-parity-nbl-ci` claim, just on a
different acceleration axis.

## What this claim does NOT say

- It does not assert that the GPU path is currently correct. The
  50K battle test gives strong *survival* evidence (99.1% physics-
  consistent) and the Rust kernel tests give strong
  *kernel-internal* evidence, but the end-to-end Python parity
  check is the missing piece.
- It does not block release. CPU is the reference path; users
  worried about GPU correctness can disable the `cuda` feature
  and fall back to the asserted CPU path.

## How to close

1. Create `tests/test_gpu_cpu_energy_parity.py` following the
   `cross_path_parity` template, parametrised over CHARMM19+EEF1
   and AMBER96.
2. Promote this claim from `kind: reference` to `kind: measurement`
   with structured tolerances populated from the cited test.
3. Update `last_verified` on first replay with the observed max
   GPU↔CPU divergence.
