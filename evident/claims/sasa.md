# SASA vs Biopython

Operational case writeup for the SASA claims in `claims/sasa.yaml`. Audience:
a proteon developer or reviewer reading the manifest who needs the why
behind the tolerances and oracle choice.

## Problem

Shrake-Rupley SASA is a well-defined integral, but each implementation makes
small choices — atomic-radii table, default probe radius, sphere
discretisation — that shift the answer by a fraction of a percent. The trust
question is therefore not "does it work?" but "how far does proteon drift
from a community-accepted reference, and is that drift bounded across the
input distribution we actually serve?"

Downstream features (RSA-driven core/surface labels, supervision-export
features, hydrophobic-burial energetics) inherit whatever bias proteon's
SASA carries, so a percent-level drift here is not cosmetic.

## Trust Strategy

Validation only. There is no proof; understanding goes only as far as
"Shrake-Rupley with consistent radii and probe should agree to a few
percent." The discipline lives in the oracle comparison.

- **Oracle**: Biopython's `Bio.PDB.SASA.ShrakeRupley`. Independently
  implemented, widely cited, shares no code with proteon.
- **Tier split**: a CI claim pins a tight tolerance on a single fixed input
  (1crn) so PRs catch regressions in seconds; a release claim pins a
  distributional tolerance on a 1000-PDB sample so the deployment shape is
  exercised before tagging.

## Evidence

- **CI**: `tests/test_sasa.py::TestBiopythonOracle::test_total_sasa_matches_biopython`
  asserts <2% relative difference on crambin total SASA. Runs on every
  push; Biopython is `pip install biopython`-cheap so no skip gating.
- **Release**: `validation/run_validation.py --n-structures 1000` records
  per-structure relative deltas in `validation/results.json`. Most recent
  reported median: 0.17%. Pass/warn/fail counts on the full corpus are
  961/2/37 (37 fails are loading-stage cascades, not SASA bugs).

Both layers exercise the same proteon entry points (`proteon.total_sasa`,
`proteon.atom_sasa`) so the CI claim is a strict subset of the release
claim's input space.

## Assumptions

- Probe radius and sphere-point class match between proteon and Biopython
  at the chosen call (`probe=1.4 A`, defaults).
- The Biopython version pinned in the proteon dev environment is the same
  one that produced the cited 0.17% number; an upstream change in
  Biopython's atomic-radii table would shift the oracle.
- 1crn passing at 2% does not generalise to all chain sizes — the release
  claim exists because of this.

## Failure Modes

- **Single-oracle blind spot.** Only Biopython is wired in. FreeSASA is in
  the install table (`tests/oracle/README.md`) but no test currently calls
  it. A shared Shrake-Rupley convention bug would be invisible to this
  case as it stands. Closing the gap is a small contribution and would
  let us move from "Oracle Comparison" toward "Cross-Tool Agreement".
- **Corpus drift.** `validation/pdbs/` is not version-pinned inside the
  manifest. Regenerating the corpus without rerunning the validation
  script would leave the cited 0.17% stale while
  `validation/results.json` still exists at the same path.
- **Forensic gap.** `results.json` records summary counts and per-structure
  relative diffs, not absolute SASA values. Replaying a divergence to root
  cause requires the corpus PDBs preserved alongside the JSON.

## Lessons

- Two tiers are worth the extra YAML. The CI claim catches regressions;
  the release claim defends the deployment shape. Collapsing them to one
  would either gate PRs on a 1000-structure run or ship without
  distributional evidence.
- Tolerances belong in the manifest, not in folklore. The 2% / 0.5%
  numbers are deliberate — they encode discretisation noise plus a margin.
  Anything larger should provoke investigation, not assertion-loosening.
- Single-oracle oracle comparisons are honest only when the gap is named.
  FreeSASA is in the failure-modes list specifically so a reviewer reading
  the claim cannot mistake "agrees with Biopython" for "agrees with the
  consensus Shrake-Rupley answer."
