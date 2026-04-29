# Fold preservation: proteon AMBER96 vs GROMACS AMBER96 L-BFGS

Operational case writeup for the AMBER96 fold-preservation release
claim in `claims/fold_preservation_amber_gromacs.yaml`.

## Problem

The CHARMM-side fold-preservation claim
(`claims/fold_preservation_charmm.md`) compares against OpenMM
CHARMM36+OBC2 — different force fields by construction. The AMBER96
side has a stronger comparator available: **GROMACS** also runs
AMBER96 with an L-BFGS energy minimiser. So the AMBER96
fold-preservation claim can be sharper — same FF, different
implementation, on the same corpus.

The trust question: across a 1000-PDB sample, does proteon's AMBER96
minimisation preserve fold geometry comparably to GROMACS's AMBER96
+ L-BFGS minimisation, at matched cutoffs and minimiser budgets?

## Trust Strategy

Validation. Same FF on both sides means this is closer to a parity
claim than the CHARMM/CHARMM36 case, but minimiser-implementation
gaps (proteon's L-BFGS vs GROMACS's L-BFGS) plus pre-processing
differences (PDBFixer for proteon, `gmx pdb2gmx` for GROMACS) keep
it on the comparator side rather than strict parity.

Pipelines:
- proteon: `validation/tm_fold_preservation_amber.py` (the AMBER
  counterpart of the CHARMM script).
- GROMACS: `validation/tm_fold_preservation_gromacs.py`. Per-
  structure: `gmx pdb2gmx -ff amber96 -water none -ignh` →
  `editconf -box 30 30 30` → `grompp` (steep EM, rcoulomb=rvdw=14
  nm) → `mdrun -nt 1` → extract CA, compute
  `proteon.tm_score(ca_pre, ca_post)`.

GROMACS settings tuned to match OpenMM/proteon as closely as
possible: `integrator = l-bfgs`, `nsteps = 100`, `emtol = 10.0`,
`rcoulomb = 14.0`, `rvdw = 14.0`, `cutoff-scheme = Verlet`,
`constraints = none`.

## Evidence

For each of 1000 sampled PDBs (or however many of the
`gmx_fold_preservation/pdbs_1k` corpus exist, see assumptions),
both pipelines:
1. Build the AMBER96 topology (PDBFixer or `pdb2gmx`).
2. Extract CA pre-minimisation.
3. Run L-BFGS minimisation, 100 steps, matched cutoffs.
4. Extract CA post-minimisation.
5. Compute `proteon.tm_score(ca_pre, ca_post)`.

GROMACS-side artifact:
`/scratch/TMAlign/proteon/validation/gmx_fold_preservation/tm_fold_gromacs.jsonl`.
Per-structure record: `pdb`, `tm_score`, `rmsd`, `n_ca`, `wall_s`,
optional `error` with parsed `Fatal error` from `gmx` stderr.

Comparison is statistical (median, p95) like the CHARMM claim. The
AMBER/AMBER setup tightens the expected gap considerably; large
divergence between proteon and GROMACS medians on the same FF is a
real problem signal, not "different FFs by design".

## Assumptions

- The corpus at
  `/scratch/TMAlign/proteon/validation/gmx_fold_preservation/pdbs_1k`
  is the input to both sides. Some PDBs may fail GROMACS prep
  (`gmx pdb2gmx` is stricter than PDBFixer about residue templates,
  HETATMs, terminal hydrogens); the GROMACS run skips those and
  records the failure in JSONL. Comparison medians are computed
  over the GROMACS-success subset.
- The `GMX` env var points at a working GROMACS binary
  (default `/scratch/TMAlign/gromacs-2026.1/build/bin/gmx`).
- proteon and GROMACS both use AMBER96 parameter values from the
  same upstream lineage; small parameter-file vintage differences
  (1-4 vdW scaling, charge dictionary version) are absorbed by the
  comparator tolerance.
- Vacuum minimisation on both sides; no implicit solvent (`-water
  none`, no OBC).

## Failure Modes

- **GROMACS prep failures bias the corpus.** `gmx pdb2gmx` rejects
  PDBs that PDBFixer accepts. The remaining "GROMACS-success
  subset" may be cleaner / smaller / less diverse than the full
  1000. The GROMACS JSONL records failures, but comparison
  medians are computed over an implicitly filtered set.
- **Different prep pipelines on the two sides.** PDBFixer (proteon)
  and `gmx pdb2gmx` (GROMACS) handle missing atoms, terminus
  hydrogens, and tautomers differently. A regression in either
  prep tool would shift one side without affecting the other.
- **Different minimisers despite same name.** Both sides nominally
  run L-BFGS; the underlying convergence criteria, line-search,
  and Wolfe conditions can differ. Matching `emtol` and `nsteps`
  bounds the budget but not the trajectory.
- **GROMACS box-size choice.** `editconf -box 30 30 30` puts every
  structure in a 30 nm cubic box; proteon has no equivalent
  notion. For very large structures (>30 nm extent) this would
  fail, but those are filtered upstream of the corpus.
- **Corpus is local-only.** Unlike the CHARMM claim corpus
  (`/globalscratch/`), the GROMACS corpus is at
  `/scratch/TMAlign/proteon/validation/gmx_fold_preservation/pdbs_1k`
  inside the repo path tree but not in the source tree. Same
  pinning issue — corpus_sha placeholder until tarballed.
- **GROMACS oracle currency.** `gromacs-2026.1` is pinned in the
  default `GMX` path but the manifest's `pinned_versions[GROMACS]`
  is a placeholder.

## Lessons

- Same-FF different-implementation is a stronger comparator than
  different-FF SOTA. The CHARMM claim absorbs the FF gap; this
  claim removes it and sharpens the test.
- Prep-pipeline diversity is itself a comparator. PDBFixer and
  `gmx pdb2gmx` accepting different subsets of PDBs is data — the
  intersection of accepted structures is a more honest "common
  ground corpus" than either pipeline's superset.
