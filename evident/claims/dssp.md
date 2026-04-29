# DSSP vs pydssp

Operational case writeup for the DSSP CI claim in `claims/dssp.yaml`.

## Problem

Proteon ships a Rust port of the Kabsch-Sander 1983 DSSP algorithm
(secondary structure assignment from CA distances + H-bond energies).
Downstream consumers (supervision export, fold preservation, search
features) inherit whatever drift proteon's DSSP carries on per-residue
assignments and on gross helix/strand fraction.

The trust question is "does proteon's DSSP agree with an independent
implementation of the same paper, on a corpus that exercises both
mostly-helix and mostly-sheet folds, well enough that downstream
features are not biased by an undetected DSSP regression?"

## Trust Strategy

Validation. Two independent implementations of one paper are the right
oracle shape: a divergence is either a proteon bug or a convention
gap, and the cause is almost always readable from the failure diff.

- **Oracle**: pydssp (Shintaro Minami, AIST) — an independent
  NumPy/PyTorch reimplementation of Kabsch-Sander. Different code,
  different research group, same paper. `pip install pydssp`.
- **Convention gap**: pydssp emits a 3-class alphabet (H / E / -),
  proteon emits the full 8-class DSSP alphabet (H, G, I, E, B, T, S,
  C). The test collapses proteon to 3-class before comparing
  (H/G/I→H, E/B→E, T/S/C→-). The collapse is **lossy on helix flavor
  and isolated-vs-ladder beta** but exactly spans pydssp's output —
  any surviving disagreement is a real H-bond detection or boundary
  difference, not an alphabet artifact.
- **Fixture corpus**: 5 structures spanning mixed (1crn, 1ubq), pure
  helix bundle (1enh), large alpha/beta (1ake), and 4-chain
  hemoglobin (4hhb). Sizes 46–574 residues.

## Evidence

`tests/oracle/test_dssp_oracle.py::TestDsspOracle` runs four
parametrised checks:

- `test_length_matches`: residue counts must agree exactly. A
  divergence here is a fixture problem (HETATM/chain filtering),
  not a scoring tolerance, so it fails hard.
- `test_3class_agreement`: per-residue exact equality on the
  collapsed 3-class alphabet, asserted at ≥95%. Empirically 97-100%
  on this corpus; the residual sits on helix/strand boundary residues
  where the H-bond energy crosses the −0.5 kcal/mol threshold and
  tiny coordinate differences flip the call.
- `test_helix_fraction_parity`: gross |H_proteon − H_pydssp| < 0.10.
  Catches systematic H-bond drift that 3-class agreement might mask
  on structures with very few helices.
- `test_strand_fraction_parity`: symmetric, |E_proteon − E_pydssp| <
  0.10.

## Assumptions

- pydssp's installed version is the same one that produced the cited
  agreement numbers; an upstream change in pydssp's H-bond detection
  thresholds would shift the oracle.
- Both tools see the same residues — proteon and pydssp must filter
  HETATMs and chains identically. The `test_length_matches` hard
  assertion guards this.
- The 95% threshold absorbs boundary-residue coin flips, not
  systematic regressions. Empirical agreement runs 97-100%; falling
  below 96% should provoke investigation, not assertion-loosening.
- The 8→3 class collapse is exact: every 8-class symbol maps to a
  unique 3-class symbol with no tie-breaking. Disagreements after
  collapse therefore reflect upstream H-bond differences only.

## Failure Modes

- **Helix-flavor and isolated-bridge resolution lost.** Collapsing
  H/G/I → H means a proteon regression that misclassifies alpha vs
  3_10 vs pi-helix is invisible. Similarly E/B → E hides isolated
  beta bridges merging with extended strand. The 8-class oracle
  (`mkdssp` C++ binary) listed in `tests/oracle/README.md` as a
  "candidate oracle" would close that gap; it is not yet wired in.
- **Single-oracle blind spot.** Only pydssp is exercised. A shared
  H-bond convention bug between proteon and pydssp would be invisible
  to this claim. The candidate `mkdssp` would also help here.
- **Corpus drift unpinned.** STRUCTURES is a Python literal, not
  pinned with a SHA in the manifest. Adding or replacing a structure
  changes the asserted set silently.
- **95% is a per-structure threshold.** A claim that passes on 5
  structures at exactly 95% has weaker evidence than one passing at
  99%. The manifest does not record the *observed* agreement number;
  `last_verified.value` will once the replay runner is wired.
