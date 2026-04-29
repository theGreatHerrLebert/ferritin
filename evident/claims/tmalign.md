# TM-align vs C++ USAlign

Operational case writeup for the TM-align CI claim in
`claims/tmalign.yaml`. Audience: a proteon developer or reviewer
deciding whether the proteon TM-align port is faithful to the canonical
C++ reference.

## Problem

Proteon ships a Rust port of TM-align (`proteon-align/src/core/`) that
the broader pipeline (alignment, search, supervision export, fold-
preservation benchmarks) inherits. The trust question is *not* "does
TM-align run?" — it is "does proteon's TM-align report TM-score, RMSD,
and aligned-residue count that match the canonical C++ implementation
on inputs the user is likely to throw at it, within bounds that reflect
real numerical-portability noise rather than masking a port bug?"

C++ TMalign and USAlign are compiled with `-ffast-math`, which reorders
floating-point operations. Rust forbids that reordering. So we cannot
expect bit-identical agreement; the right tolerance is "agreement to
4-5 decimal places on TM-score, ±0.05 Å on RMSD, ±2 residues on
alignment length" — which is what the tests assert.

## Trust Strategy

Validation. The implementation understanding is decent (proteon's port
is documented in `proteon-align/src/core/`, parameter conventions and
the `d0` formula match), but the discipline that makes the port
defensible is the per-pair comparison against `USalign` on a curated
fixture set.

- **Oracle**: the C++ USAlign binary, located via `$USALIGN_BIN` or the
  `USalign` name on `$PATH`. The test skips silently when neither is
  resolvable, so dev loops without USAlign installed remain green;
  proteon CI is expected to provide the binary.
- **Fixture corpus**: 7 PDB pairs covering self-alignment (`1ubq vs
  1ubq`), low-similarity (`1ubq vs 1crn`), and same-fold pairs (`1lyz
  vs 2lzm`, `1tit vs 2gb1`, `5pti vs 1bpi`, `1mbn vs 2hbg`). Each pair
  is a `pytest` parametrised case; available pairs are filtered by
  file-existence at collection time.

## Evidence

`tests/oracle/test_tmscore_oracle.py::TestTMScoreVsCpp` runs C++
USAlign in `-outfmt 2` mode, parses the tabular output, and asserts:

- `rust.tm_score_chain1 == approx(cpp.tm2, abs=5e-4)` and the
  symmetric chain-2 case. **Note the convention swap** — proteon's
  `tm_score_chain1` is normalised by chain-2 length, which matches
  USAlign's `TM2`. The test handles this; downstream consumers that
  index by "chain1" should be aware of it.
- `rust.rmsd == approx(cpp.rmsd, abs=0.05)`.
- `abs(rust.n_aligned - cpp.n_aligned) <= 2`.

A separate `TestFastMode::test_fast_vs_normal` checks that proteon's
fast mode (`fast=True`) lands within 80–110% of the non-fast TM-score.
That is **internal** parity (proteon vs proteon), not oracle parity,
and is intentionally not part of this claim — see failure modes.

## Assumptions

- The USAlign binary is byte-stable across CI invocations; we have not
  pinned a specific commit in the manifest and rely on the host system
  to provide a recent build.
- The fixture PDBs at `test-pdbs/` (some inside the proteon repo, some
  at the sibling `/scratch/TMAlign/test-pdbs/`) do not change without
  re-running this test. We do not pin `corpus_sha` because the corpus
  spans two locations and the proteon repo SHA pins only one.
- Proteon's coordinate parsing agrees with USAlign's at the input
  layer; that's covered by a separate I/O parity claim. Disagreement
  there would surface as TM-score drift here.

## Failure Modes

- **Corpus drift unpinned.** A regression that adds, removes, or
  replaces a fixture would silently expand or shrink the asserted set
  without updating the claim. A future tightening would either pin a
  manifest of fixture SHAs or move the fixtures wholly into the proteon
  source tree and use the source SHA as the pin.
- **Fast-mode out of scope.** `TestFastMode` is a within-proteon parity
  check; it should be modeled as a separate cross-path-parity claim
  alongside other fast/slow path comparisons. Combining it here would
  conflate "agrees with USAlign" with "agrees with itself" and dilute
  the claim's queryability.
- **Independent-implementation gap.** USAlign and proteon's TM-align
  share algorithmic ancestry; both descend from the original C++
  TMalign. Strictly speaking, USAlign is an *independent
  implementation* of the same algorithm, not an *independent oracle*.
  A stronger Rule 5 ("Use Independent Validation") posture would add a
  second comparator (the original TMalign C++ binary, or a Python
  re-implementation) so that a shared numerical ancestor cannot mask a
  shared bug.
- **Oracle skip on missing binary.** The test silently skips when
  USAlign is not on `$PATH`. CI must positively assert that the
  binary is present, or a stack with a missing oracle could pass the
  test suite while exercising none of this claim.
