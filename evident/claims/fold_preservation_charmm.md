# CHARMM19+EEF1 fold preservation vs OpenMM CHARMM36+OBC2

Operational case writeup for the release-tier
`proteon-charmm19-fold-preservation-vs-openmm-release-1k-pdbs` claim.

## Problem

Per-component force-field parity (`charmm19-vs-ball-corpus-50k`) tells
you the **gradients agree** with a reference implementation.
Fold-preservation tells you the **integrated gradient over a
minimization trajectory agrees**. The two are different: a small
per-component diff can compound over 100 steps into a meaningful fold
shift, and a large per-component diff can cancel out into negligible
fold motion. Both signals are necessary; neither is sufficient.

For CHARMM19+EEF1 specifically, BALL is the canonical per-component
oracle, and OpenMM CHARMM36+OBC2 is the natural fold-preservation
oracle — different parameter set, different solvation model, but the
same modeling target ("CHARMM-family force field that preserves protein
folds during energy minimization on small wall-clock budgets").

## Approach

1. **Same input corpus** on both arms: a 1000-PDB random sample (seed=42)
   from `validation/pdbs/`. Each side runs PDBFixer preprocessing
   independently.
2. **Same minimizer budget**: 100 L-BFGS-style steps. Neither side runs
   to convergence — the question is "does the structure survive
   minimization without unfolding?", not "did we find the minimum?".
3. **Same TM-score function**: proteon's TM-align port computes the
   per-PDB TM-score on **both** arms, so the comparison is on the
   physics, not on alignment idiosyncrasies.
4. **Per-PDB pairing**: results joined by `pdb` into a single canonical
   JSONL (`charmm_pair_1k.jsonl`) with both sides' TM-score, RMSD,
   energies, and a pre-computed `tm_diff = openmm − proteon`.

## Outcome

Headline (dataset 2026-04-12/13):

| Side | median TM | mean TM | RMSD median |
|---|---|---|---|
| proteon CHARMM19+EEF1 | **0.9944** | 0.9875 | 0.45 Å |
| OpenMM CHARMM36+OBC2 | 0.9991 | 0.9950 | 0.21 Å |
| `openmm − proteon` TM diff | median **+0.0040** | | |

n_ok = 886 / 1000 attempted. 14 % per-side failures (proteon-side
free-CYS HG topology, OpenMM-side non-finite first-step forces on
PDBFixer-reconstructed atoms). The 86 % pass rate is the comparison
surface.

Both medians sit above 0.99 — fold is preserved on the typical
structure. The +0.004 diff favours OpenMM, which is expected: the
OBC2 implicit solvent produces a smoother gradient field than EEF1's
context-aware solvation, so OpenMM "settles" into a lower-RMSD
neighbourhood on the minimization budget. The diff is well inside
the claim's 0.01 tolerance band.

## Why this isn't a "CHARMM19 implementation correctness" claim

CHARMM19 ≠ CHARMM36, and EEF1 ≠ OBC2. The diff captures **both** the
parameter set difference (CHARMM19 vs CHARMM36) and the implicit
solvent difference (EEF1 vs OBC2). It is **not** "proteon's CHARMM19
implementation is wrong by 0.004 TM-score" — it is "two valid
implementations of fold-preserving CHARMM-family minimization agree
to 0.004 TM-score at the median".

Per-component CHARMM19 implementation correctness is the corpus oracle
claim's job (`charmm19-vs-ball-corpus-50k`); fold preservation is what
this claim covers.

## What this claim isn't

- **Not a benchmark of speed.** proteon's CHARMM19+EEF1 minimization
  was ~30 × faster than OpenMM CHARMM36+OBC2 on the same hardware,
  but speed is not in the tolerance set here. A separate
  performance claim could surface that.
- **Not a fold-prediction claim.** Fold preservation tests that
  energy minimization doesn't *destroy* an existing fold. Predicting
  folds from sequence is a different problem entirely.
- **Not gated on every CI run.** This is release-tier; it runs as
  part of release-bundle locking on monster3, not on every PR. The
  CI-tier gate is the per-component CHARMM-vs-BALL crambin check.
