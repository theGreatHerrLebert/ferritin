# 50K production-scale battle test

Operational case writeup for the 50K-PDB battle-test release claim
in `claims/scale_50k_battle_test.yaml`.

## Problem

The per-component oracle claims (SASA, DSSP, AMBER96, OBC GB,
TM-align, …) all fix the *correctness* of individual subsystems on
small fixtures or 1000-PDB samples. None of them answer a simpler,
more existential question users actually ask:

> *Will proteon survive my 50,000-PDB pipeline run, and produce
> physically reasonable output across the distribution of crystal
> structures it might see in the wild?*

This claim is about **survival at scale + physics-consistency**, not
about parity with an external implementation. It is the largest
single piece of operational evidence proteon has, and not having it
in the EVIDENT manifest leaves the trust story incomplete on
exactly the dimension users care about most.

## Trust Strategy

Validation against a **physical invariant** at scale, not against an
external oracle. The "oracle" in this claim is physics itself: a
folded protein under a reasonable force field should minimise to a
**negative** total energy. Any structure where minimisation finishes
with `final_energy >= 0` is either pathological input (unresolvable
clashes, broken topology) or a bug in the pipeline.

Two complementary metrics define "survival":
- **Coverage**: fraction of input PDBs that finish with
  `status=ok` (no exception, no skip). At 50K random PDBs from the
  PDB this is non-trivial — broken templates, unusual residues,
  malformed records, etc. all eat structures.
- **Physics consistency**: among `status=ok` structures, fraction
  with `final_energy < 0`. This is the proxy for "minimisation
  actually did the right thing" — energies should go negative as
  van-der-Waals + electrostatic + bonded relax to the local minimum.

## Evidence

`validation/stage3_50k_gpu_results.jsonl` — 47,183 records, one per
PDB, with schema `{pdb, atoms, status, skipped, initial_energy,
final_energy, steps, converged, exception?}`.

`validation/stage3_50k_gpu.log` — driver console log capturing
chunk-by-chunk progress, final summary, total wall time.

**Run parameters** (extracted from the log header):
- 47,183 PDBs from `/globalscratch/dateschn/proteon-benchmark/pdbs_50k`
- chunks of 500 PDBs
- `minimize_steps=50` (half the Stage-1 smoke default of 100,
  optimised for throughput at scale)
- GPU: NVIDIA RTX 5090, CUDA Compute 12.0, auto-dispatch via
  proteon's `cuda` feature
- 128 CPU cores for the parallel host-side orchestration

**Outcome** (from log tail summary):
- Total: 47,183 PDBs
- Errors: 1,697 (3.6%) — mostly `load_error` for malformed inputs
- Skipped: 244 — pre-loading filters (size, content)
- Processed (`status=ok`): 45,242 (95.9% of input)
- Negative final energy: 44,837 / 45,242 = **99.1% physics-consistent**
- Wall time: 12,631 s = **3.5 h on RTX 5090 + 128 cores**
- Throughput: 3.6 structures/s

The Stage-1 smoke test (`benchmark/stage1_smoke.py`) asserts
`pct_negative_total >= 90%` on a 100-PDB sample as the green/red
gate. This claim re-asserts that invariant at 470× the scale and
on the GPU code path.

## Assumptions

- The `/globalscratch/dateschn/proteon-benchmark/pdbs_50k` pool is
  the corpus the cited 99.1% number was produced against. Pool
  drift between runs would invalidate the cited rate without
  signal; `corpus_sha` is `PENDING-PIN-50K-PDB-POOL-TARBALL` until
  tarballed.
- 47,183 ≠ 50,000 because the pool is named "50k" but contains
  fewer files (filtering at acquisition time). The manifest's `n`
  reflects the actual pool size, not the marketing name.
- `final_energy < 0` is a meaningful proxy for "minimisation did
  the right thing." For folded proteins under CHARMM19+EEF1 at
  reasonable parameters, this holds robustly on the corpus
  distribution. It does **not** generalise to unfolded peptides,
  designed proteins, or unusual oligomeric states.
- 50 minimise steps is enough to reach a negative energy on most
  reasonable inputs. The Stage-1 smoke test uses 100 steps as
  default; 50 steps at 50K scale is a deliberate
  performance-vs-coverage tradeoff.
- The GPU code path (CHARMM19+EEF1 energy + SASA on CUDA) is what
  was exercised. CPU fallback paths exist but were not the subject
  of this run; CPU correctness is covered by the per-component
  oracle claims.

## Failure Modes

- **Driver script is not in the source tree.** The artifact is
  preserved (`validation/stage3_50k_gpu_results.jsonl`,
  `validation/stage3_50k_gpu.log`) but the actual driver that
  produced them is *not currently committed*. The structure
  matches a Stage-3 sibling of `benchmark/stage1_smoke.py`
  (same JSONL schema, same negative-energy invariant, scaled with
  parallelism + GPU), but a future replay requires re-creating the
  driver from the parameters in the log header before the next
  release. The `command` in this manifest is the **reconstructed**
  invocation, not a literal script path that exists today.
- **Survival != correctness.** A pipeline that produces *plausible
  but wrong* energies on every structure could pass this claim.
  The per-component oracle claims (SASA vs Biopython, AMBER96 vs
  OpenMM, OBC GB vs OpenMM) cover correctness; this claim covers
  survival. Both are needed.
- **Physics invariant has tail cases.** Small peptides, designed
  proteins, or structures with severe steric clash that
  minimisation cannot resolve in 50 steps may legitimately end at
  positive energy. The 99.1% rate already absorbs this, but a
  regression that pushes structurally fine PDBs into the positive-
  energy bucket would manifest as a drop below 95%.
- **Single-hardware claim.** Cited on RTX 5090. A regression
  specific to older CUDA compute capabilities (sm_70, sm_80) or
  to CPU fallback would not surface here. Memory cites
  GPU↔CPU 1e-11 parity on crambin but no asserting test wires it
  at the public API; that gap belongs in a separate cross-path
  parity claim.
- **Duration is real.** 3.5 h on RTX 5090 + 128 cores is the
  literal cost of running this claim. Release-tier replay needs
  a budgeted slot.
- **Skip + error counts not gated.** The 95.9% `status=ok` rate is
  reported but not asserted; only the negative-energy rate among
  `status=ok` is gated. A regression that increased the `ok` →
  `error` ratio (e.g., a loader regression on edge-case PDBs)
  would show in the log but not break the structured tolerance.

## Lessons

- Survival-at-scale and per-component-correctness are **different
  claims**, not different evidence for the same claim. Conflating
  them weakens both. The 50K test does not say "AMBER96 is
  correct"; it says "the pipeline does not crash on 50K random
  PDBs." Both statements are useful, separately.
- The driver-not-in-source gap is itself a process lesson: scripts
  that produce load-bearing artifacts need to live in the source
  tree at the same SHA as the artifact, or the artifact's
  reproducibility is permanently compromised. This is the single
  most important hardening this claim demands of the EVIDENT
  workflow.
