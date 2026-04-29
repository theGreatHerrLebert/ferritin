# SASA cross-tool vs FreeSASA (gap)

Operational case writeup for the FreeSASA cross-tool *reference*
claim in `claims/sasa_freesasa.yaml`.

## Why this is a gap, not a measurement

The existing SASA claims (`proteon-sasa-vs-biopython-ci` and
`proteon-sasa-vs-biopython-release-1k-pdbs`) are single-oracle:
Biopython's `Bio.PDB.SASA.ShrakeRupley` is the only reference. Per
EVIDENT Rule 5 ("Use Independent Validation: agreement alone is
not sufficient"), single-oracle claims have a real failure mode —
a shared Shrake-Rupley convention bug between proteon and
Biopython would be invisible.

`tests/oracle/README.md` lists FreeSASA in the install table. It
is `pip install freesasa`-cheap, no source build required. But:
**it is not currently installed in the proteon dev venv** (verified
2026-04-29) and **no test calls it**. The gap is fixable in an
afternoon; the cost of *not* fixing it is the failure-modes
section of every existing SASA claim citing the gap explicitly.

## Required test design

A new oracle test, either extending `tests/test_sasa.py` with a
`TestFreeSASAOracle` class or creating
`tests/oracle/test_sasa_freesasa_oracle.py`, that:

1. Imports FreeSASA via `pytest.importorskip("freesasa")`.
2. Computes total SASA on the 1crn fixture using both proteon's
   default and FreeSASA's default settings.
3. Asserts |proteon - freesasa| / freesasa < 2% (mirroring the
   Biopython tolerance).
4. Optionally extends to per-atom SASA agreement at a slightly
   looser tolerance (FreeSASA exposes per-atom output cleanly).

The release-tier 1000-PDB benchmark (`validation/run_validation.py`)
should also gain a FreeSASA branch alongside the existing Biopython
one, recording per-structure deltas in `validation/results.json`.

## Trust strategy

Validation by **cross-tool agreement** (the next pattern up from
"Oracle Comparison" in `evident/patterns/README.md`). Three
independent Shrake-Rupley implementations (proteon, Biopython,
FreeSASA) agreeing within tolerance is materially stronger evidence
than two — the failure mode "all three implementations share the
same wrong convention" is implausible enough to discount.

## What this claim does NOT say

- It does not invalidate the existing Biopython SASA claims; those
  remain correct, just incomplete by Rule 5.
- It is not an upgrade in tolerance — both Biopython and FreeSASA
  use the same Shrake-Rupley algorithm with similar atomic-radii
  conventions. The 2% tolerance is the right band on both axes.

## How to close

1. `pip install freesasa` in the proteon `.venv`.
2. Add a `TestFreeSASAOracle` class to `tests/test_sasa.py` (or a
   new oracle file).
3. Add the FreeSASA branch to `validation/run_validation.py::test_sasa`.
4. Promote this claim from `kind: reference` to `kind: measurement`
   with structured tolerances.
5. Update the existing two SASA claims' failure_modes to remove the
   "FreeSASA gap" note (or downgrade it to "FreeSASA covered in
   sibling claim").
