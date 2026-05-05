# Fold preservation at 50K: proteon CHARMM19+EEF1 vs OpenMM CHARMM36+OBC2

Operational case writeup for the 50K-scale fold-preservation
release-tier claim in `claims/fold_preservation_charmm_50k.yaml`.
Layered on top of the 1k claim (`claims/fold_preservation_charmm.yaml`):
the 1k claim establishes that the two minimizers preserve fold to within
0.005 TM-score at the median across a small random sample; this claim
establishes that the agreement holds across two orders of magnitude more
structures.

## Why this claim exists at 50K

The 1k version's headline (proteon median TM=0.9944, OpenMM median
TM=0.9991, tm_diff = +0.0040) is well below the 0.01 band, but the 1k
sample size is too small to distinguish "consistent agreement" from
"sample-noise that happens to land inside the band". The 50K extension
gives the claim statistical power: a real implementation drift would
show as a tail-shift across many structures even if it stays under the
median band on any 1k sample.

This claim also serves as the production-default force field's
cross-implementation oracle. proteon's CHARMM19+EEF1 is the production
default for fold preservation, structure preparation, and downstream
analysis; agreement at 50K with OpenMM's CHARMM36+OBC2 (different
parameter set, but the same minimization outcome under the same scoring
metric) is the strongest cross-tool evidence proteon ships for the
production force field's correctness.

## Approach

1. Same 50 000-PDB random sample of wwPDB used by other 50K claims.
2. Both arms run a 100-step L-BFGS-style minimizer at effective
   NoCutoff; neither runs to convergence (the metric is fold survival,
   not minimum-finding).
3. proteon-side: `validation/tm_fold_preservation.py` with `N_PDBS=50000`.
4. OpenMM-side: `validation/tm_fold_preservation_openmm.py` with
   `N_PDBS=50000`. PDBFixer skip-missing-atoms applied (PR #54).
5. Joiner: `validation/fold_preservation/join_fold_preservation.py`
   reads from `PROTEON_OUTPUT_DIR` and emits the canonical
   `charmm_pair_1k.jsonl` artifact (the filename keeps the `1k` suffix
   for renderer compatibility but the records reflect the actual N).
6. TM-score on both arms uses proteon's TM-align port — alignment
   noise cancels in the diff.

## Outcome

PENDING — first 50K run scheduled. Headline numbers will be filled in
at first lock.

Expected (extrapolating from the 1k headline):

```
n_attempted ~ 50 000
n_ok        ~ 40 000-45 000  (10-20% population narrowing for missing
                              heavy atoms + non-protein residues)
proteon median TM   ~ 0.99
OpenMM  median TM   ~ 0.99
tm_diff median      ~ +0.004
```

Substantial deviation from these would indicate either a population
bias (the 50K sample skews different from the 1k sample) or a real
implementation drift that emerges at scale.

## What this claim isn't

- **Not a force-field correctness claim.** The two arms use *different*
  force fields (CHARMM19+EEF1 vs CHARMM36+OBC2). Agreement on TM-score
  after minimization is an "is the minimum well-defined" claim, not a
  "do the energies agree" claim.
- **Not a per-structure regression test.** Per-structure CI is in
  `tests/test_fold_preservation_invariants.py`.
- **Not a CHARMM force-field implementation oracle.** That role is
  `claims/forcefield_charmm19_ball_corpus_50k.yaml` (proteon vs BALL,
  same CHARMM19+EEF1 parameter set on both arms).
