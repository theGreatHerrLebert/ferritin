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
   claim. Wall time: ~225 minutes (v0.1.4) on a 32-worker
   AMD Threadripper-class host. Wall time grew vs v0.1.3 because
   pebble runs each task in its own subprocess (no cascade
   short-circuit).

## Outcome

Headline (v0.1.4, lock time 2026-05-04):

```
n_attempted=44 210
n_ok        = 5 309
n_skip      = 37 258  (18 912 missing-heavy-atoms,
                       18 344 BALL NaN/inf typer,
                       2 non-AA residues)
n_err       =  1 643  (1 604 60s timeouts,
                       35 PDB-parse ValueError,
                       4 BALL SIGSEGV — isolated by pebble)

                  median   p95    p99    pass/5309   band
bond_stretch       0.06%  2.21%  8.80%   4731/5309   <1.0%
angle_bend         0.01%  0.39%  2.82%   5179/5309   <1.0%
vdw                0.01%  0.13%  1.47%   5268/5309   <2.5%
electrostatic      0.02%  0.02%  0.02%   5289/5309   <1.0%
solvation          1.61%  2.16%  2.39%   5304/5309   <5.0%
torsion            0.00%  0.00%  0.00%   5309/5309   <2.5%
improper_torsion   0.70% 23.11% 200.34%       —      <2.5% (median only)
```

All seven components meet their per-tier band on the median.
improper_torsion's median (0.70%) is well within band; the tail
(p95 ~23%, p99 ~200%) is the residue-template variant coverage
gap surfaced at scale, the same gap visible in the 1k claim.
bond_stretch's p99 ~8.8% is a similar tail (free-CYS HG and
disulfide-CYS edge cases).

## Why the n_ok jumped 807 → 5 309 between v0.1.3 and v0.1.4

The v0.1.3 bundle hit `concurrent.futures.ProcessPoolExecutor`'s
"broken pool" cascade: a single worker SIGSEGV from BALL on a
pathological topology propagated to all 31 in-flight futures, even
with `max_tasks_per_child=1`. The cascade rate climbed to ~93% on
a single pass against the 50K random wwPDB sample.

PR #44 migrated the runner to `pebble.ProcessPool.schedule(timeout=…)`,
which runs each task in its own subprocess and handles SIGSEGV as
a per-task failure with no pool-level damage. The 4 isolated
SIGSEGVs in the v0.1.4 run are now per-task records; without
pebble they would have re-cascaded.

The remaining 37 258 skips are now dominated by well-defined,
scientifically documented out-of-population categories rather than
runner-quality cascades — see PR #47 for the missing-heavy-atoms
skip rationale (PDBFixer.addMissingAtoms hangs deterministically on
certain inputs).

## What this claim isn't

- **Not a per-PDB regression test.** The crambin gate
  (`claims/forcefield_charmm19_ball.yaml`) covers per-PDB drift.
- **Not a fold-preservation claim.** That's
  `claims/fold_preservation_charmm.yaml` — different metric (TM
  score after minimisation), different oracle (OpenMM CHARMM36
  + OBC2).
- **Not a "100% PDB coverage" claim.** It's a 5 309-records-out-of-
  44 210 claim against a population narrowed to "well-resolved
  wwPDB structures BALL can type", framed honestly. Widening the
  population requires upstream fixes in BALL's typer and PDBFixer's
  reconstruction loop, not force-field changes.
