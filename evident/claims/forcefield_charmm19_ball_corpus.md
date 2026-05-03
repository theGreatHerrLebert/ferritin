# CHARMM19+EEF1 vs BALL on a 1000-PDB random sample

Operational case writeup for the 1k corpus CHARMM19-vs-BALL release-tier
claim in `claims/forcefield_charmm19_ball_corpus.yaml`. Layered on top
of the per-component crambin gate
(`claims/forcefield_charmm19_ball.yaml`) — the unit oracle establishes
parity on a single fixture; this claim establishes that the parity
generalises across a random sample of the wwPDB.

## Problem

The crambin oracle gates CHARMM19+EEF1 at the per-component level on
a single 46-residue fixture. That's necessary but not sufficient: the
crambin composition (1 PHE, 2 TYR, 6 disulfide CYS, no histidines, no
free CYS) doesn't exercise the residue-template variants, terminal
edge cases, or large-loop geometries the production-default force
field has to handle.

A 1000-PDB random sample from `validation/pdbs/` (a curated proteon
corpus seeded for benchmark reproducibility) is the smallest corpus
that exposes the residue-template coverage gaps in a statistically
useful way: PHE/TYR ring exclusions, disulfide CYS impropers, free
CYS HG handling, terminal NH3+/COO- variants, and the long tail of
multi-chain biological assemblies are all sampled.

## Approach

1. **Filter**: every PDB in the 1k sample is checked by PDBFixer; any
   structure that fails to clean (missing residues that can't be
   modelled, nucleic acids the protein-only filter let through) is
   recorded as `skipped` in the JSONL and excluded from headline
   statistics.
2. **Prepare**: PDBFixer removes hetero, replaces nonstandard residues,
   and adds missing heavy atoms; proteon then writes the polar-H PDB.
   BALL reads the same polar-H PDB with `add_hydrogens=False` so both
   tools evaluate the same atoms — H placement is no longer a source
   of disagreement.
3. **Compute**: proteon's `compute_energy(ff="charmm19_eef1",
   nonbonded_cutoff=1e6)` and `ball.charmm_energy(pdb, use_eef1=True,
   nonbonded_cutoff=1e6)` run NoCutoff. The adapter extracts BALL's
   pure electrostatic via `nonbonded - vdw - solvation` (the older
   `nonbonded - vdw` form silently double-counted EEF1; see PR #24).
4. **Compare**: per-component relative diff
   `|proteon - BALL| / |BALL|` for each of the 7 components, with
   per-tier tolerance bands matching the crambin unit claim.

## Outcome

5 of 7 components passing on the headline run (commit `c6a1c5d`,
2026-04-30): bond_stretch, angle_bend, vdw, electrostatic, and
solvation all meet the per-tier band on the median. proper_torsion
and improper_torsion were x-failed at the time of authorship — both
have since closed in PR #26 and are now passing on the
crambin gate; this claim's `last_verified.value` will reflect that
once a fresh corpus run lands.

## What this claim isn't

- **Not a single-PDB regression test.** The crambin gate
  (`claims/forcefield_charmm19_ball.yaml`) is the per-PDB regression
  surface. This claim runs on a sample; small per-PDB drift doesn't
  fail it.
- **Not a 50K-scale claim.** That's
  `claims/forcefield_charmm19_ball_corpus_50k.yaml`, which surfaces
  the residue-template coverage tail at scale.
- **Not a pebble-isolated run.** ProcessPoolExecutor with
  `max_tasks_per_child=1` is the runner; pathological PDBs (~14% of
  the 1k corpus) trigger BrokenProcessPool cascades that are
  documented in `failure_modes`. The pebble migration is tracked in
  the 50K claim's failure_modes.
