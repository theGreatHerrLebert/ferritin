# MSA vs MMseqs2 (search + result2msa pipeline)

Operational case writeup for the MSA research-tier claim in
`claims/msa.yaml`.

## Problem

Proteon ships a clean-room port of the MMseqs2 search path
(prefilter → ungapped → gapped Smith-Waterman → MSA assembly) under
`proteon-search/`. The port targets AlphaFold-style MSA feature
generation: same query, same target DB, same parameters → "the same
MSA, modulo documented convention gaps."

The trust question is delicate: MMseqs2 is the upstream the port is
*derived from*, so byte-for-byte agreement on the MSA file is too
strict (would mask any legitimate engineering improvement). At the
same time, "the MSAs look similar" is too loose for an ML pipeline
where each row is a feature vector. The right shape is a **structured
shared-hit comparison** with explicit tolerance bands and a named
convention gap.

## Trust Strategy

Validation by oracle on the operational unit (the MSA), not on the
implementation (CIGAR-by-CIGAR equivalence). Two comparison modes:

- **Positional** — compare MSA rows by rank (both sorted by alignment
  score descending). Sensitivity differences between proteon's
  prefilter and upstream's mean most rows compare different targets;
  this metric is reported for completeness but not gated.
- **Shared-hit** — match MSA rows by target DB key, compare only
  targets found by both pipelines. Isolates the alignment/assembly
  math from prefilter divergence. **This is the gated claim.**

Per the validation script's documented results (50 queries vs
uniref50_1m, 2026-04-19): shared-hit residue agreement median 0.896,
deletion agreement median 0.974, gap agreement median 0.710 (not
gated). All 50 queries PASSED the contract.

## Evidence

`validation/msa_oracle.py`:

1. Run proteon search + MSA on each query against the target DB.
2. Run upstream `mmseqs search` + `mmseqs result2msa` on the same
   query + DB with matched parameters.
3. For each query: find the shared-hit set (target DB keys present
   in both MSAs), compare row-by-row residue / deletion / gap
   agreement, write per-query metrics to
   `validation/msa_oracle_report.json`.
4. Aggregate medians, means, p5 across queries.

## The convention gap

Diagnosed 2026-04-19 on A0A0C1M9X2 (1262 residues, 171 shared hits
against uniref50_1m): **0/171 shared hits had identical CIGARs.**

Root cause: proteon's Smith-Waterman extends alignments further into
low-identity tails. Typical case: proteon aligns positions 6–1262 vs
upstream 298–658. Query-start is close (median 4 residues apart);
query-end diverges substantially (62/171 differ by >20 residues).

Diagnosis: alignment-policy difference, not math error.

- Both tools implement standard SW with affine gap penalties.
- Proteon uses more aggressive extension on low-identity regions.
- All affected hits are remote homologs (20–39% fident).
- Within the overlapping aligned region, residues match ~90% and
  deletions match ~97%.

This is the same pattern as the AMBER96 NoCutoff convention gap
documented in `devdocs/ORACLE.md` §"When the oracle is also wrong,
case 3": a policy choice that creates a measurable gap, not a math
bug. The shared-hit residue/deletion agreement metrics are calibrated
to be sensitive to the math (residue identity, deletion projection)
without breaking on the boundary policy.

## Assumptions

- The query fixture (`proteon-search/tests/data/oracle_fixture.fasta`)
  and target DB (`uniref50_1m`) are byte-stable; corpus_sha is a
  PENDING-PIN placeholder until tarballed.
- Upstream MMseqs2 (`/globalscratch/aweissen/bio/MMseqs2/build/bin/mmseqs`)
  is the version that produced the cited 2026-04-19 numbers; an
  upstream sensitivity-default change would shift the prefilter
  hit-set without affecting the shared-hit metrics.
- The shared-hit comparison subset is large enough per query (median
  171 hits on the cited diagnostic) to bound the per-query agreement
  metric tightly.
- Sensitivity-independence is verified within noise (s=5.7 vs s=7.5
  produced identical shared-hit metrics on the same corpus).

## Failure Modes

- **Tier=research, not gated in CI.** This claim runs on demand
  during dev cycles and before tags. A regression that ships in CI
  but breaks the MSA oracle would not be caught until the next
  manual run. Promotion to release tier needs (a) the corpus SHA
  pinned, and (b) the runtime budgeted for a release pass.
- **Hit-set Jaccard expected low.** Proteon's prefilter is much more
  sensitive than upstream at default settings (10–50× more hits at
  `-s 5.7`/`-s 7.5`). Hit-set Jaccard sits at 0.01–0.10 by design;
  the shared-hit comparison is what the claim gates on. A reviewer
  might mistake the low Jaccard for a regression.
- **MSA depth ratio not gated.** Proteon saturates `max_seqs` while
  upstream returns 2–256 hits. The depth ratio is reported but
  asserts nothing.
- **Convention gap could mask drift.** The 0.75 residue-agreement
  threshold accommodates the boundary-extension gap; a real drift
  in the alignment math that produces ~75% residue agreement would
  pass. Tighter thresholds would require closing the boundary
  policy gap.
- **Single 50-query batch.** Cited results are on 50 queries; the
  oracle script can scale but the manifest's `n` reflects only the
  default smoke run.

## Lessons

- Operational-unit oracle (MSA file) > implementation-detail oracle
  (CIGAR equivalence). The latter would have failed on every shared
  hit with no actionable signal; the former gives clean structured
  metrics that point at boundary policy when they fail.
- Document the convention gap in the case writeup so future
  reviewers reading "0.75 residue agreement is a low-looking
  threshold" understand it is a deliberate policy concession, not
  an aspirational target.
