# DSSP 8-class vs mkdssp (gap)

Operational case writeup for the 8-class DSSP *reference* claim in
`claims/dssp_8class.yaml`.

## Why this is a gap, not a measurement

The existing DSSP claim (`proteon-dssp-vs-pydssp-ci`) collapses
proteon's 8-class output (H, G, I, E, B, T, S, C) to pydssp's 3-class
(H, E, -) before comparison. The collapse rules:

- H, G, I → H (all helix flavors)
- E, B → E (extended strand + isolated bridge)
- T, S, C → -

Proteon emits 8-class data in production. Downstream consumers
that index on helix flavor (alpha vs 3₁₀ vs π) or on isolated-bridge
vs ladder beta receive **uncovered** output: the manifest currently
makes no claim about the correctness of those distinctions.

`tests/oracle/README.md` explicitly lists `mkdssp` (the C++
reference implementation, Kabsch & Sander's own descendant) as a
candidate oracle that would close this gap. It is not currently
wired in.

## Required test design

A new oracle test, `tests/oracle/test_dssp_mkdssp_oracle.py`, that:

1. Detects `mkdssp` via `MKDSSP_BIN` env var (same pattern as
   `REDUCE_BIN`); skip cleanly if unset.
2. Runs `mkdssp` on the same 5-structure corpus as the pydssp
   oracle (1crn, 1ubq, 1enh, 1ake, 4hhb).
3. Parses `mkdssp` output (DSSP plain-text format) into the 8-class
   alphabet.
4. Asserts per-residue agreement with proteon's 8-class output at
   ≥95% (boundary residues will still disagree, but helix-flavor
   and isolated-bridge calls should mostly match).

Two sub-claims emerge naturally:

- **Helix-flavor distribution** — fraction of H + G + I cells, with
  per-flavor agreement (H vs G vs I).
- **Isolated bridge** — fraction of B (vs E) cells, agreement on
  the bridge/ladder split.

Both are output fields proteon emits today but the manifest does
not currently claim correctness on.

## Trust strategy

Validation by independent implementation. `mkdssp` is the C++
reference maintained by Maarten Hekkelman (rooted in the original
Kabsch-Sander codebase). It is not pip/conda-installable on most
modern toolchains; the test would either build it from source in
CI or skip when the binary is unavailable.

This is the same install-difficulty class as `reduce` — the
EVIDENT claim accepts the env-gate convention.

## What this claim does NOT say

- It does not assert that the 3-class collapse (the existing pydssp
  claim) is wrong. That claim remains correct for users who only
  need helix/strand/loop labels.
- It does not block release. The 3-class claim covers the dominant
  use case; this gap matters specifically for users consuming the
  full 8-class output.

## How to close

1. Build `mkdssp` from source under `/scratch/TMAlign/dssp/` (sibling
   pattern with `reduce`).
2. Create `tests/oracle/test_dssp_mkdssp_oracle.py`.
3. Add `MKDSSP_BIN` to `tests/oracle/README.md` install table.
4. Promote this claim to `kind: measurement` with structured
   per-flavor tolerances.
