# AMBER96 vs OpenMM, 1000-PDB single-point oracle

Operational case writeup for the AMBER96-vs-OpenMM release claim in
`claims/forcefield_amber_openmm.yaml`.

## Problem

OpenMM is proteon's **authoritative** AMBER96 oracle. The BALL CI
claim pins component-level math on a single fixture; this claim
answers a different question: *across a representative 1000-PDB
sample, does proteon's AMBER96 single-point total energy agree with
OpenMM's at a tight per-structure relative tolerance, or are there
structural classes where the two diverge?*

Without the distributional check, a regression that shifts crambin by
~0.5% (passes BALL's 1% band) but blows up on long flexible loops
would ship invisibly.

## Trust Strategy

Validation. OpenMM is the de-facto reference AMBER96 implementation
in the molecular-dynamics community: maintained by Stanford's Pande
group lineage, used by AlphaFold, OpenForceField, and most modern
ML/MD pipelines. Two independent tools (proteon, OpenMM) loading the
**same prepared structure** (PDBFixer-prepped to identical atoms +
hydrogens) and reporting the same single-point energy is the
strongest practical AMBER96 oracle.

- **Oracle**: OpenMM (`ForceField("amber96.xml")`,
  `NoCutoff`, no constraints, no rigid water). CPU platform with
  Threads=1 for determinism.
- **Corpus**: 1000 PDBs sampled (seed=42) from the
  `proteon-benchmark/pdbs_50k` 50K random PDB pool. Same pool used
  by the fold-preservation benchmarks.
- **Prep**: PDBFixer for both tools — find/fix missing residues and
  atoms, add hydrogens at pH 7.0, strip waters and heterogens. Both
  tools then load the same prepared PDB so the only remaining
  variable is each tool's energy implementation.

## Evidence

`validation/amber96_oracle.py`:
- For each of 1000 sampled PDBs, runs PDBFixer prep, OpenMM single-
  point at NoCutoff, proteon single-point with `nbl_threshold=10**9`
  to force the exact O(N²) path.
- Records per-structure absolute and relative energy difference,
  per-force-group OpenMM breakdown, per-component proteon
  breakdown, and walltime.
- Writes JSONL artifact at
  `/globalscratch/dateschn/proteon-benchmark/amber96_oracle.jsonl`.
- Prints summary stats: `|diff|` mean / median / p95 / max in
  kJ/mol, and `rel` mean / median / p95 / max.

The script's docstring states the success criterion as
`per-structure |E_proteon - E_openmm| / |E_openmm| < 1e-3 on total`.
That is the contract; the artifact lets a reviewer verify per-
structure compliance.

## Assumptions

- PDBFixer's prep is deterministic enough that proteon and OpenMM
  see byte-identical atoms and positions on each structure;
  re-running PDBFixer with the same input produces the same output
  (modulo float-stable RNG-free passes).
- The 50K PDB pool at `/globalscratch/dateschn/proteon-benchmark/`
  is the corpus the cited result was produced against; corpus_sha
  in the manifest is a `PENDING-PIN-*` placeholder until the pool
  is canonicalised (tarball hash or git-LFS pin).
- OpenMM's `amber96.xml` ships the parameter file proteon's AMBER96
  implementation also targets; a unit-conversion or charge-table
  difference would systematically skew every structure.
- `NoCutoff` on both sides means the comparison is independent of
  proteon's default 15 Å cutoff policy. The release claim does not
  speak to performance-tier energies; that is a separate claim about
  the policy gap.

## Failure Modes

- **Per-structure success rate not encoded.** The script's docstring
  states the per-structure threshold but the runner does not gate
  pass/fail per structure — it dumps everything to JSONL and a
  reviewer reads the medians. Until the manifest's `last_verified`
  starts recording the observed >=95% pass rate, a reader cannot
  distinguish "median 0.2%, all <1e-3" from "median 0.2%, 5% of
  structures at 5%". The structured tolerance below pins the median;
  the per-structure rate sits in prose only.
- **Corpus is external and unpinned.** The 50K PDB pool lives on a
  scratch filesystem outside the repo. No SHA, no manifest, no
  tarball. A regenerated pool with subtly different chain filtering
  would produce different numbers without any signal.
- **PDBFixer is part of the trust chain.** Hydrogen placement, missing-
  atom rebuild, and heterogen stripping all happen in PDBFixer
  before either tool sees the input. A PDBFixer regression would
  shift both proteon and OpenMM in lockstep — invisible to this
  claim.
- **Component-level partitioning differs between tools.** OpenMM
  reports per-`ForceGroup`; proteon reports per-component. Bond,
  angle, vdW, coulomb, improper map cleanly between the two but
  the script asserts only on the total. A compensating-component
  regression would be visible in the JSONL but not gated.
- **NoCutoff only.** Proteon's production default uses a 15 Å
  cutoff for performance; the gap between NoCutoff and 15 Å is a
  documented policy choice and would surface as a separate
  performance-vs-accuracy claim.

## Lessons

- "Authoritative" oracle != "automatically authoritative result".
  PDBFixer in the prep chain is a dependency that shifts the
  goalposts; a future tightening should pin the PDBFixer version in
  `pinned_versions` (or replace it with proteon's own prep
  pipeline) to reduce that surface.
- Distributional release tier is not a substitute for
  per-component CI. The BALL claim and this OpenMM claim measure
  different things; both are needed.
