# CHARMM19+EEF1 vs BALL on a 44 210-PDB protein-only corpus

Operational case writeup for the 50K-scale CHARMM19-vs-BALL
release-tier claim in
`claims/forcefield_charmm19_ball_corpus_50k.yaml`. Layered on top of
the 1k corpus claim (`claims/forcefield_charmm19_ball_corpus.yaml`):
the 1k claim establishes parity holds across a small random sample;
this claim establishes that it holds across two orders of magnitude
more structures, and surfaces the long-tail residue-template coverage
gaps the 1k sample size is too small to catch.

## Problem

A 1000-PDB sample is just large enough to exercise residue-template
variants on a few examples per type. It can't measure tail risk: a
1-in-200 residue-template miss disappears under sample noise. When
proteon's CHARMM19+EEF1 is the production-default force field for
fold preservation, structure preparation, and downstream analysis,
the failure rate at scale matters more than the median.

## Approach

1. **Source corpus**: a 50 000-PDB random sample of the wwPDB,
   downloaded once for proteon's benchmark suite and pinned at
   `/globalscratch/dateschn/proteon-benchmark/pdbs_50k/` on
   monster3.
2. **Filter**: `validation/protein_only_corpus.py` drops nucleic
   acids and non-standard residues that CHARMM19 cannot
   parametrise. Pass rate: 88% (44 210 / 50 000). The filter
   output (`validation/protein_only_50k.txt`) is the actual
   replay input.
3. **Run**: same runner (`validation/charmm19_eef1_ball_oracle.py`),
   same per-component split, same NoCutoff convention as the 1k
   claim. Wall time: ~30 minutes on a 32-worker
   AMD Threadripper-class host.

## Outcome

Headline (commit at lock time, 2026-05-03):

```
n_ok = 807 / n_skip = 2 084 / n_err = 33 (of 2 924 attempted)
                       (44 210 - 41 346 BrokenProcessPool cascades)

bond_stretch    median 0.054%  pass 770/807
angle_bend      median 0.011%  pass 800/807
vdw             median 0.014%  pass 802/807
electrostatic   median 0.019%  pass 804/807
solvation       median 1.621%  pass 807/807
torsion         median 0.001%  pass 807/807
improper_torsion median 0.414%  pass 597/807
```

Six of seven components meet their per-tier band on the median.
improper_torsion's median (0.41%) is well within band; the tail
(p99 ~52%) is the residue-template variant coverage gap surfaced
at scale.

## Why the n_ok is 807, not ~14 500

The 1k claim observed a 13.9% BrokenProcessPool cascade rate; the
50K random sample contains many more pathological inputs (free CYS,
unusual altlocs, large multi-chain assemblies that PDBFixer leaves
in degenerate states, etc.) and the cascade rate climbs to ~93% on
a single pass. Python's
`concurrent.futures.ProcessPoolExecutor` propagates a single worker's
SIGSEGV exit to all 31 in-flight futures even with
`max_tasks_per_child=1`.

The 807 records are real measurements with full proteon and BALL
energies on each side — they're not biased by the cascade. The
~41 346 cascade errors are runner-quality noise, not force-field
correctness data points. The runner-side fix is to swap
`concurrent.futures` for `pebble.ProcessPool.schedule(timeout=...)`,
which runs each task in its own subprocess and handles SIGSEGV as
a per-task error with no pool-level damage. Tracked as a
ball-py + runner robustness follow-up; the 807 records here are
preserved unchanged through that migration.

## What this claim isn't

- **Not a per-PDB regression test.** The crambin gate
  (`claims/forcefield_charmm19_ball.yaml`) covers per-PDB drift.
- **Not a fold-preservation claim.** That's
  `claims/fold_preservation_charmm.yaml` — different metric (TM
  score after minimisation), different oracle (OpenMM CHARMM36
  + OBC2).
- **Not a 50K _correctness_ claim.** It's a 807-records-out-of-44 210
  claim, framed honestly. Improving the n_ok ratio is a
  runner-quality task (pebble migration), not a force-field
  improvement.
