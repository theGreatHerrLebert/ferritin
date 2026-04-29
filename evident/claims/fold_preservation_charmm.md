# Fold preservation: proteon CHARMM19+EEF1 vs OpenMM CHARMM36+OBC2

Operational case writeup for the fold-preservation release claim in
`claims/fold_preservation_charmm.yaml`.

## Problem

Proteon's value proposition for structural-bioinformatics workflows
includes "minimisation that doesn't wreck the fold." Single-point
energy parity (covered by the BALL and OpenMM oracle claims) is
necessary but not sufficient — minimisation drives a structure
through gradient steps that compound any per-step bias. A
force-field that produces correct energies but drives crambin's CA
trace 5 Å away from the input is not usable for ML pipelines that
expect minimised input to remain a faithful representation of the
deposited fold.

The trust question: across a representative 1000-PDB sample, does
proteon's CHARMM19+EEF1 minimisation preserve fold geometry — as
measured by per-structure pre-vs-post-minimisation TM-score — at a
level comparable to OpenMM's CHARMM36+OBC2 SOTA reference, run on
the same input PDBs under matched minimiser budgets?

## Trust Strategy

Validation by **comparator**, not strict parity. Proteon and OpenMM
use **different force fields** (CHARMM19+EEF1 vs CHARMM36+OBC2);
per-PDB TM-scores are not expected to match exactly. The claim is
distributional: proteon's median fold-preservation TM should be
high in absolute terms (>= 0.99) and within a small margin of the
SOTA OpenMM number on the same corpus.

Two parallel benchmarks:
- `validation/tm_fold_preservation.py` — proteon CHARMM19+EEF1
  minimisation, writes `tm_fold_preservation.jsonl`.
- `validation/tm_fold_preservation_openmm.py` — OpenMM
  CHARMM36+OBC2 + LocalEnergyMinimizer, writes
  `tm_fold_preservation_openmm.jsonl`.

Both runs use the same 1000-PDB sample (seed=42 from the 50K pool),
PDBFixer prep on both sides for matched topology, and matched
minimiser budgets (proteon `minimize_steps=100, tol=0.1 kcal/mol/Å`;
OpenMM `LocalEnergyMinimizer.minimize(tolerance=10 kJ/mol/nm,
maxIterations=100)`).

Per the documented validation status: proteon median TM = 0.9945,
OpenMM CHARMM36+OBC2 median TM = 0.9991. Proteon is ~30× faster
wall-time and ~1.5× more tolerant of raw PDBs.

## Evidence

For each of 1000 sampled PDBs, both pipelines:
1. PDBFixer prep (find/fix missing residues, add H at pH 7).
2. Extract CA coordinates pre-minimisation.
3. Build the FF system and run the minimiser.
4. Extract CA coordinates post-minimisation.
5. Compute TM-score(pre, post) using `proteon.tm_score`.

Output JSONL columns: `pdb`, `n_ca_pre`, `n_ca_post`,
`initial_energy_kj`, `final_energy_kj`, `tm_score`, `rmsd`,
`n_aligned`, `wall_s`, optional `error`.

The artifact is the *pair* of JSONLs: comparison is statistical
(median, p95, distribution shape), not per-row.

## Assumptions

- The 1000-PDB sample (seed=42) is the same corpus on both sides.
  Any drift in `PDB_DIR` between runs would silently de-pair the
  comparison.
- PDBFixer prep is byte-identical for both pipelines per PDB; both
  scripts call `PDBFixer.findMissingAtoms / addMissingAtoms /
  addMissingHydrogens(7.0)`.
- OpenMM CHARMM36+OBC2 is the right SOTA comparator for fold
  preservation in implicit solvent. CHARMM19+EEF1 (proteon) is the
  united-atom equivalent; SOTA-vs-classic gap is expected to be
  small but non-zero.
- 1000 structures is enough to bound the median tightly. p95 and
  tail behaviour are visible in the JSONL but not gated by the
  structured tolerance.
- TM-score on CA pre-vs-post is the right fold-preservation metric.
  RMSD is co-recorded but not used in the structured claim — large
  RMSD with high TM happens routinely on flexible loops.

## Failure Modes

- **Different force fields in the comparator.** This is not a
  parity claim; OpenMM is a *contextual SOTA reference*, not the
  oracle saying "proteon's CHARMM19 result is correct." A
  reviewer mistaking it for parity would over-interpret the
  evidence.
- **Two separate runs, no strict pairing.** The median comparison
  works because both runs cover the same PDB set, but per-PDB
  comparison is not asserted. A regression that shifts a few PDBs
  by a lot while leaving the median intact would not surface here.
- **Corpus is external and unpinned.**
  `/globalscratch/dateschn/proteon-benchmark/pdbs_50k` is on a
  scratch filesystem; corpus_sha is `PENDING-PIN`. Until pinned,
  the cited median is reproducible only if the pool itself is
  preserved.
- **Minimiser-budget mismatch.** Both sides run 100 steps with
  matched tolerances, but the underlying minimisers differ
  (proteon's L-BFGS vs OpenMM's `LocalEnergyMinimizer` which is
  also L-BFGS). A subtle convergence-criterion difference could
  produce systematic gap that looks like a force-field difference.
- **PDBFixer is in the trust chain on both sides.** A PDBFixer
  regression would shift both medians in lockstep — invisible
  here. Pinning the PDBFixer version in `pinned_versions` is the
  right tightening.
- **Cutoff differs.** OpenMM uses `CutoffNonPeriodic` at 1.0 nm
  with implicit OBC2; proteon CHARMM19+EEF1 has its own implicit-
  solvent cutoff policy. The comparison absorbs this as part of
  the "different FFs" assumption rather than asserting a matched
  cutoff.

## Lessons

- Comparator claims earn their keep when an absolute SOTA reference
  is unattainable. CHARMM19+EEF1 cannot bit-match CHARMM36+OBC2 by
  construction, but OpenMM-on-the-same-corpus is the right way to
  argue "proteon is in the same league."
- The 30× speed gap is real and worth recording but does not belong
  in the same claim. A separate research-tier claim or benchmark
  should pin throughput; this claim is about correctness.
