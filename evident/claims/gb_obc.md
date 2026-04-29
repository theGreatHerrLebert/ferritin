# OBC GB vs OpenMM on crambin

Operational case writeup for the OBC-GB-vs-OpenMM release claim in
`claims/gb_obc.yaml`.

## Problem

Proteon's OBC (Onufriev-Bashford-Case) generalised-Born implicit
solvent landed in 2026-04-15 (Phase B). It implements the GB pair
sum with OBC1 parameters (ε_in=1, ε_out=78.5, offset=0.09 Å, α=0.8,
β=0, γ=2.909125) and integrates with the AMBER96 vacuum energy to
produce the AMBER96+OBC total.

The trust question for this component is sharp: *does proteon's
single-point AMBER96+OBC total agree with OpenMM's, and does the GB
component itself match OpenMM's GB single-point, on a controlled
fixture (crambin) with identical atoms and hydrogens fed to both
tools?*

## Trust Strategy

Validation. OpenMM's `amber96_obc.xml` (which adds `GBSAOBCForce` to
`amber96.xml`) is the consensus reference for AMBER+OBC implicit
solvent in the MD community. The script forces matched conditions
across both tools so the only remaining variable is each tool's
implementation:

- **Force field**: AMBER96 on both sides.
- **Implicit solvent**: OBC1 with explicit parameter pin in the
  module docstring.
- **Cutoff**: `NoCutoff` on both, with proteon's
  `nbl_threshold=10**9` forcing the exact O(N²) path.
- **Hydrogens**: PDBFixer is run once and the prepped PDB is fed to
  both proteon and OpenMM, eliminating H-placement as a degree of
  freedom.

## Evidence

`validation/amber96_obc_oracle.py`:

1. PDBFixer-prep crambin to a temp PDB.
2. OpenMM single-point with `amber96.xml + amber96_obc.xml` →
   `total_kj`. Decompose by removing the GB force from a parallel
   system → `vacuum_kj`. Subtract → `gb_kj`.
3. Proteon single-point with `ff="amber96_obc"` → `total_kj`,
   `vacuum_kj = total - solvation`, `gb_kj = solvation`.
4. Compare: `delta_total`, `rel_total`, `delta_gb`, `rel_gb`.

The script's documented contract (the "Phase D PASS" condition):
`rel_total < 1e-2 (1%) AND rel_gb < 5e-2 (5%)`. Phase A baseline
(GB stub returning zero) was the earlier failure state; Phase B
landing closes the contract.

Per memory, the achieved numbers on crambin: ≤5% GB / ≤1% total vs
OpenMM, GPU matches CPU to 1e-11.

## Assumptions

- PDBFixer prep is deterministic so proteon and OpenMM see
  byte-identical atoms and hydrogen positions.
- Proteon's `ff="amber96_obc"` selects the OBC1 parameter set
  matching OpenMM's `amber96_obc.xml`. A different OBC variant
  (OBC2, GBn, GBn2) would produce a 5-15% drift even with correct
  math.
- The `amber96.xml + amber96_obc.xml` pair in OpenMM is the
  combination producing the cited reference. Loading
  `amber96_obc.xml` alone (without the bonded/non-bonded base) is
  not equivalent and would fail at templating.
- The 1% total / 5% GB tolerances are calibrated for crambin's
  size and dielectric profile. Larger or more solvent-exposed
  proteins may require recalibration; the manifest scopes the
  claim explicitly to crambin.

## Failure Modes

- **Single fixture.** Crambin only, n=1. A regression on larger
  proteins or proteins with unusual SASA distributions would not
  surface. Extending to a small distributional release-tier
  benchmark (e.g. 100 PDBs) is a reasonable tightening once the
  contract holds at scale.
- **OBC1 only.** Proteon may grow OBC2/GBn variants later; this
  claim does not speak to those.
- **Released-not-tightened tolerances.** The 5% GB band is wider
  than ideal; per memory, achieved is ≤5%. An OpenMM update or a
  proteon Born-radius regression that creeps inside the band
  without crossing it would not be caught. The manifest's
  `last_verified.value` should record the observed `rel_gb` so
  drift below the band becomes visible.
- **Phase A vs Phase B/C/D legacy.** The script's `if rel_total <
  1e-2 and rel_gb < 5e-2` block returns 0 (pass) on Phase B+ but
  also returns 0 (pass) on Phase A "stub" because it short-circuits
  on `gb_kj < 1e-6`. Now that Phase B has landed, that
  short-circuit should be removed; until it is, the script could
  silently pass on a regressed-to-zero GB. CI must check that
  `gb_kj` is non-zero for the claim to be honest.

## Lessons

- Pin the parameter regime in the case writeup, not just in the
  test code. OBC has multiple variants; OBC1 is what proteon
  targets and the manifest needs to make that explicit so a future
  switch to OBC2 surfaces as a manifest change rather than a
  tolerance creep.
- Phase contracts age. `validation/amber96_obc_oracle.py` was
  written as the Phase D contract while Phase A was live; now that
  Phase B has landed, the Phase A short-circuit is a foot-gun and
  should be removed. The EVIDENT manifest is the right place to
  flag that obligation.
