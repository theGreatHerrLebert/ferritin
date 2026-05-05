# Fold preservation at 50K: proteon AMBER96 vs OpenMM AMBER96

Operational case writeup for the 50K-scale fold-preservation AMBER
release-tier claim in `claims/fold_preservation_amber_50k.yaml`.
Sister claim to the CHARMM 50K version
(`claims/fold_preservation_charmm_50k.md`); same approach, different
force field, but with one critical difference: **both arms use the same
parameter set (AMBER96)**, so any disagreement is purely implementation,
not parameter-set divergence.

## Why this claim exists at 50K

This is the cleanest cross-implementation force-field oracle proteon
ships at 50K scale. The 1k version's tm_diff median of +0.0028 (proteon
TM=0.9959, OpenMM TM=0.9992) is the tightest cross-tool agreement
across all proteon's release-tier oracles, exactly because the two
arms aren't being asked to agree on a force-field parameter set —
they're being asked to agree on how to integrate the same parameter
set toward a local minimum.

A 50K extension turns this from "the implementations agree on a small
random sample" into "the implementations agree across the wwPDB
population". For AMBER96 specifically, that's the strongest possible
evidence proteon's AMBER96 minimization is correct: if a 50K random
sample shows median tm_diff staying inside the 0.005 band, the
implementations have to be doing the same physics on the same input
positions.

## Approach

Identical to the CHARMM 50K version (see
`fold_preservation_charmm_50k.md`), with the substitutions:

- proteon-side: `validation/tm_fold_preservation_amber.py` (proteon
  with `ff="amber96"`, `constrain_heavy=False`)
- OpenMM-side: `validation/tm_fold_preservation_openmm_amber.py`
  (OpenMM with `amber96_obc.xml`)

Both runners now honour `N_PDBS` env so a 50K invocation is a one-line
override of the existing 1k command. Skip-missing-atoms (PR #54)
applies on the OpenMM side; the proteon side handles its own per-PDB
errors via `batch_prepare`.

## Outcome

PENDING — first 50K run scheduled.

Expected (extrapolating from the 1k headline):

```
n_attempted ~ 50 000
n_ok        ~ 40 000-45 000
proteon median TM   ~ 0.996
OpenMM  median TM   ~ 0.999
tm_diff median      ~ +0.003
```

If the 50K headline lands in this range, AMBER96 minimization in
proteon is implementation-equivalent to OpenMM AMBER96 across the
production-relevant population. That's the trust signal v0.2.0 is
built on.

## What this claim isn't

- **Not a CHARMM correctness claim.** AMBER96 is not proteon's
  production-default force field (CHARMM19+EEF1 is).
  `claims/forcefield_charmm19_ball_corpus_50k.yaml` and
  `claims/fold_preservation_charmm_50k.yaml` carry that role.
- **Not a single-point energy claim.** That's
  `claims/forcefield_amber_openmm_50k.yaml` — different metric,
  different bug class (per-term coefficient mismatches that don't
  show up after minimization).
