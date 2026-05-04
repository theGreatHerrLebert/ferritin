# AMBER96 fold preservation: proteon vs OpenMM

Operational case writeup for the release-tier
`proteon-amber96-fold-preservation-vs-openmm-release-1k-pdbs` claim.

## Why this is the cleanest force-field oracle proteon has

Both arms claim **AMBER96** specifically — same parameter set, same
nominal physics. The only thing different is the implementation:
proteon's Rust + PyO3 connector vs OpenMM's reference C++. So the
per-PDB TM-score diff captures **pure implementation drift**, not
parameter-set divergence.

Compare with the CHARMM fold-preservation claim
(`fold_preservation_charmm`), which compares CHARMM19+EEF1 vs
CHARMM36+OBC2 — different parameter set AND different solvent
treatment. That diff captures both the physical model difference
and the implementation. It's still a useful claim but a noisier
oracle than this one.

## Approach

Same shape as the CHARMM claim:

1. **Same input corpus**: 1000-PDB random sample (seed=42) from
   `validation/pdbs/`. PDBFixer preprocessing on each side.
2. **Same minimizer budget**: 100 L-BFGS-style steps. Fold-survival
   metric, not minimum-finding.
3. **Same TM-score function**: proteon's TM-align port computes both
   sides' TM. Diff isolates the physics.
4. **Per-PDB pairing** by `pdb` into `amber_pair_1k.jsonl` with both
   sides' values + a pre-computed `tm_diff = openmm − proteon`.

## Outcome

Headline (dataset 2026-04-13):

| Side | median TM | mean TM | RMSD median |
|---|---|---|---|
| proteon AMBER96 | **0.9959** | 0.9911 | 0.27 Å |
| OpenMM AMBER96 | 0.9992 | 0.9952 | 0.18 Å |
| `openmm − proteon` TM diff | median **+0.0028** | | |

n_ok = 864 / 1000 attempted. 13.6 % per-side failures (proteon-side
H-placement edge cases, OpenMM-side non-finite first-step forces on
PDBFixer-reconstructed atoms). 86.4 % pass rate is the comparison
surface.

The +0.003 TM diff is the **smallest cross-implementation force-field
agreement** number proteon has on a release-tier corpus. proteon and
OpenMM agree on AMBER96 fold preservation to 3 millipoints of
TM-score at the median.

## Why this isn't a CHARMM-correctness claim

AMBER96 is not proteon's production-default force field
(CHARMM19+EEF1 is). This claim is the cleanest implementation-vs-
implementation oracle proteon has, but it does NOT extend to
CHARMM19+EEF1. Production CHARMM correctness is carried by:

- **`charmm19-vs-ball-corpus-50k`** (per-component oracle vs BALL)
- **`fold_preservation_charmm`** (fold-preservation oracle vs OpenMM
  CHARMM36+OBC2)
- **`charmm19-vs-ball-ci`** (per-component crambin gate)

Together those three plus this claim form proteon's force-field trust
pyramid: per-component math (BALL) + fold preservation (OpenMM, both
CHARMM and AMBER lines) + CI regression guards.

## What this claim isn't

- **Not a benchmark of speed.** proteon's AMBER96 minimization wall
  time vs OpenMM's is not in the tolerance set here.
- **Not a fold-prediction claim.** Fold preservation tests that
  energy minimization doesn't *destroy* an existing fold; it does
  not predict folds from sequence.
- **Not gated on every CI run.** This is release-tier; runs as part
  of release-bundle locking on monster3, not on every PR.
- **Not a replay-from-image claim yet** (as of v0.1.4). The runner
  scripts have hardcoded host paths and the claim's full re-produce
  path requires the v0.2.0 runner-parametrization work tracked in
  issue #42 + #48. Verify-from-locked-JSONL works today; full
  re-produce requires that work to land.
