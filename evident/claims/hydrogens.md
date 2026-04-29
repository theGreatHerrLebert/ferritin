# Hydrogen placement vs reduce

Operational case writeup for the hydrogen-placement CI claim in
`claims/hydrogens.yaml`.

## Problem

Proteon places hydrogens during structure preparation
(`proteon-connector/src/add_hydrogens.rs`) under two regimes:
- **Polar-only** for CHARMM19 (matches CHARMM19's united-atom carbons
  that absorb non-polar C-H into inflated heavy radii).
- **Full H** for AMBER96 (every heavy atom gets its protons).

H placement upstream of force-field energies, fold preservation, and
SASA — wrong N-H or O-H positions cascade into bad H-bond detection,
bad clash counts, and slightly off SASA. The trust question is "does
proteon place the **rigid** hydrogens (backbone N-H, CA, methylene
CH₂, aromatic CH) within tens-of-milli-Ångströms of the canonical
Richardson-Lab `reduce` placement, on a small corpus that exercises
both polar-only and full-H modes?"

## Trust Strategy

Validation against the community-standard reference. `reduce`
(Richardson Lab, Duke) is the most-cited H-placement tool in the
field; its output is the de-facto canonical answer for rigid H
positions. Bonded oracle: same residue topology, different H
optimisation strategy.

- **Oracle**: `reduce` C++ binary. Gated by `REDUCE_BIN` (and
  optionally `REDUCE_DB`) env vars. Test skips silently when not set
  so dev loops without `reduce` installed remain green; CI must
  positively check.
- **Fixture corpus**: 2 structures (1crn, 1ubq) — small and well-
  characterised, deliberately scoped to the rigid H subset where the
  claim is sharp.
- **Parametrisation**: each structure runs twice, polar-only (CHARMM19
  regime) and full-H (AMBER96 regime). 2 × 2 = 4 measurement points.

## Evidence

`tests/oracle/test_reduce_hydrogen_oracle.py::TestReduceHydrogenOracle`
runs two assertions per (structure, FF mode) cell:

- `test_rigid_parity`: per-rigid-H-atom Euclidean distance ≤ 0.1 Å
  after optimal matching (Hungarian assignment within each rigid
  group). Rigid groups are defined in the test as backbone N-H, CA-H,
  methylene CH₂, and aromatic CH — positions where bond geometry
  fixes H within the heavy-atom frame.
- `test_backbone_amide_coverage`: |proteon_amide_H_count −
  reduce_amide_H_count| ≤ 2. Sanity check that the placement actually
  ran across the chain; tolerates proline/N-terminal handling
  differences without asserting them.

The test **strips hydrogens** from the input PDB before running both
tools, so neither is biased by pre-existing refined H positions.

## Assumptions

- `REDUCE_BIN` is set and points at a working `reduce_src/reduce`
  build; CI must positively assert this rather than silently skipping.
- The 0.1 Å rigid-placement tolerance is a structural floor, not a
  numerical convenience. Bond length and bond angle ambiguity at heavy
  atoms typically allows ~0.05 Å placement variance; 0.1 Å absorbs
  that with margin.
- The Hungarian-style optimal matching in `_collect_rigid_residuals`
  is the correct comparison strategy. A naive index-paired comparison
  would fail on methylene CH₂ where reduce and proteon may swap H1/H2
  labels with the same physical placement.
- 2 structures (1crn, 1ubq) are sufficient *for rigid H placement*.
  They cover backbone amides, methylenes (Lys, Arg side chains),
  aromatic CH (Phe, Tyr in 1ubq). A claim about non-rigid placement
  (rotamer-optimised polar H, methyl rotamers) would need a different
  corpus and a different tolerance.

## Failure Modes

- **Rigid subset only.** This claim says nothing about the convention
  gaps `tests/oracle/README.md` flags as "not asserted":
  - Methyl rotamers (CH₃ on Ala, Val, Leu, Ile, Thr, Met) — both tools
    place valid CH₃ but at different 3-fold orientations.
  - sp2 amide NH₂ (Asn HD21/HD22, Gln HE21/HE22, Arg NE/NH1/NH2) —
    in-plane placement convention differs by ~120°.
  - Rotatable polar H (Ser/Thr/Tyr OH, Cys SH, Lys NH3+) — reduce
    optimises by H-bond scoring, proteon places at template default.
- **Env-gated oracle.** `REDUCE_BIN` unset → test skips silently. CI
  must explicitly check that the binary is installed; a misconfigured
  pipeline can pass the test suite while exercising none of this
  claim.
- **2 structures is small.** Larger or unusual structures (membrane
  proteins, structures with non-standard residues, very high-B-factor
  loops) are not exercised. A regression specific to those would not
  surface.
- **Polar-only and full-H share most rigid H atoms.** The test
  exercises both modes but the rigid subset overlaps heavily. A
  regression that affects only non-overlapping placement (e.g. polar-
  only-mode dropping a rigid H that full-H keeps) would be partly
  visible.

## Lessons

- Define what is "rigid" precisely and assert against that, not
  against "all H atoms". The convention-gap residuals (methyl, sp2
  amide, rotatable polar) are *physically real differences*, not bugs.
  Asserting a tight tolerance on the union would force false-positive
  failures or assertion-loosening to absorb known divergence.
- Oracle tests are claim-shaped, not test-shaped. Two assertions per
  cell (rigid parity + amide coverage) rather than one big union test
  lets a regression localise to "rigid H positions" or "chain coverage"
  cleanly in CI output.
