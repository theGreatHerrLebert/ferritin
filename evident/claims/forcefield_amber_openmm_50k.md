# AMBER96 single-point energy: proteon vs OpenMM at 50K corpus

Operational case writeup for the 50K-scale AMBER96-vs-OpenMM release-tier
claim in `claims/forcefield_amber_openmm_50k.yaml`. Sister claim to the
fold-preservation AMBER pair (`claims/fold_preservation_amber.yaml` at
1k, `claims/fold_preservation_amber_50k.yaml` at 50K) — same parameter
set on both arms, but the metric is energy magnitude before minimization
rather than TM-score after minimization.

## Why this claim exists at 50K

The 1k version was retired in PR #49 in favour of the fold-preservation
pair. At 1k scale, the additional capability coverage of a single-point
energy comparison was deemed redundant when fold-preservation already
exercised both implementations end-to-end.

At 50K scale the calculus shifts: a population-level systematic
energy-magnitude drift (e.g. a per-term coefficient mismatch that
averages to zero across 1k structures but sums to a measurable bias
across 50K) is detectable in this claim but invisible in the fold-pres
TM-score metric. The two claims test different bug classes:

- **fold-pres**: detects bugs that change the *minimum* (gradient
  direction or stiffness)
- **single-point energy**: detects bugs that change the *value*
  (per-term coefficient, exclusion-list off-by-one, charge-derivation
  drift) without changing the minimum

Both are needed at the trust level v0.2.0 wants to ship.

## Approach

1. **Source corpus**: same 50 000-PDB random sample of wwPDB used by the
   CHARMM 50K oracle, filtered through `validation/protein_only_corpus.py`
   to drop nucleic acids.
2. **Runner**: `validation/amber96_oracle.py` (pebble-migrated in PR #56)
   with `PROTEON_PDB_LIST=validation/protein_only_50k.txt` and
   `N_PDBS=50000`. Per-task pebble subprocess isolation prevents the
   `BrokenProcessPool` cascade that limited the v0.1.3 CHARMM 50K run.
3. **PDBFixer prep**: skip-missing-atoms (PR #47/#54) — comparison
   surface is "well-resolved wwPDB", not "everything PDBFixer can repair".
4. **Both arms see identical atoms**: the runner writes the
   PDBFixer-prepped topology to a temp PDB and feeds it to BOTH proteon
   and OpenMM, so any disagreement is in the energy implementation, not
   the input parsing.

## Outcome

PENDING — first 50K run scheduled. Headline numbers (n_attempted, n_ok,
median rel_diff, p95/p99/max, pass rate per population narrowing
category) will be filled in at first lock.

The 1k version's headline was: median rel_diff ≈ 0.2% on total energy,
0.5% on the gated components at NoCutoff. The 50K version is expected
to track those numbers within sample-noise tolerance across the
narrower well-resolved-wwPDB population.

## What this claim isn't

- **Not a CI gate.** The per-PDB CI claim
  (`claims/forcefield_amber_ball.yaml` if restored, or the AMBER96
  invariants under `tests/test_amber_invariants.py`) covers the
  per-structure regression check.
- **Not a fold-preservation claim.** That's
  `claims/fold_preservation_amber.yaml` at 1k +
  `claims/fold_preservation_amber_50k.yaml` at 50K — different metric,
  different bug class.
