# I/O parity vs Biopython + Gemmi

Operational case writeup for the I/O parity CI claim in
`claims/io.yaml`.

## Problem

Proteon's PDB / mmCIF loader (`proteon-io/`, bridging pdbtbx) is the
input boundary for every downstream component — alignment, SASA,
energy, search, supervision export. A single-character disagreement at
this layer (chain ID, atom name, B-factor) silently propagates into
every claim that runs after it. The trust question is therefore "do
two independent reference parsers (Biopython, Gemmi) see the same
structure that proteon does, on a corpus that exercises the common
surface (single-chain PDB, mmCIF, multi-chain assembly, edge cases)?"

## Trust Strategy

Validation by **cross-tool agreement**, not just oracle comparison.
Biopython and Gemmi are independent parsers maintained by different
teams; agreeing with both narrows the convention-bug failure mode
that single-oracle SASA still has open.

- **Oracles**: Biopython (`Bio.PDB`) and Gemmi (`gemmi`), both via
  pip. Both gated by `pytest.importorskip` so dev loops without them
  remain green; CI must explicitly install both.
- **Fixture corpus**: 4 standard files in `STANDARD_FILES`
  (`tests/oracle/conftest.py`) — `1ubq.pdb`, `1ubq.cif`, `4hhb.pdb`
  (4-chain hemoglobin), `1crn.pdb` (crambin). Mixes single-chain,
  multi-chain, and PDB-vs-mmCIF.
- **Edge-case fixtures** (`models.pdb`, `insertion_codes.pdb`) are
  defined in conftest but tested separately under relaxed expectations
  — out of scope here; a future claim can pin those tolerances
  individually.

## Evidence

`tests/oracle/test_io_oracle.py` runs three classes against the
4-fixture corpus, parametrised per-fixture:

- `TestCounts`: model_count, chain_count, atom_count, residue_count,
  chain_ids — exact equality across proteon, Biopython, Gemmi.
- `TestCoordinates`: per-atom coordinates within 0.001 Å vs Gemmi
  (always) and vs Biopython (skipped when atom counts diverge due to
  alt-conformer handling differences).
- `TestMetadata`: per-atom B-factors and occupancies within 0.01 vs
  Gemmi; atom names, residue names, per-atom chain IDs exact
  equality vs Gemmi; elements with at most max(1, n//100) mismatches
  (i.e. ~1% with a floor of 1).

`TestInternalConsistency` inside the same file is a within-proteon
sanity check (numpy array views vs per-atom accessors), not an oracle
comparison; intentionally not part of this claim.

## Assumptions

- Biopython's `PDBParser`/`MMCIFParser` and Gemmi's `read_structure`
  are byte-stable across runs; no upstream regression in either
  parser.
- The 4 fixture files at `test-pdbs/` (mixed proteon-internal and
  `/scratch/TMAlign/test-pdbs/` on host) do not change without a
  regenerated test run.
- The 0.001 Å coordinate tolerance reflects the PDB format's
  3-decimal-place limit, not a numerical convention gap. mmCIF can
  carry more precision but the corpus-wide tolerance ties to the PDB
  floor.
- The element-mismatch tolerance (≤1% with a floor of 1) absorbs
  Gemmi's heavier element-guessing on records that omit the
  element column; it is **not** a license for proteon to mislabel
  elements at scale.

## Failure Modes

- **Corpus drift unpinned.** STANDARD_FILES is a 4-entry
  Python-literal list; the manifest does not pin file SHAs. Adding,
  removing, or replacing a fixture mutates the asserted set silently.
  Tightening would either move all fixtures inside the source tree
  and rely on the source SHA, or add a sidecar fixture-SHA file to
  the EVIDENT manifest's `inputs:`.
- **Two parsers, shared conventions.** Biopython and Gemmi are
  independent codebases but both implement the same wwPDB / mmCIF
  spec. A spec-level convention shift could move both oracles in
  lockstep without proteon changing.
- **Edge cases out of scope.** Multi-model NMR (`models.pdb`),
  insertion codes, and alt-conformers are defined in conftest as
  edge-case fixtures and not asserted here. A regression specific to
  those formats would not surface in this claim.
- **mmCIF coverage is one fixture.** Only `1ubq.cif` exercises the
  mmCIF parser. Asymmetric-unit-only structures and large mmCIF-only
  PDBs (10k+ atoms) are not covered.

## Lessons

- Cross-tool agreement is materially stronger than single-oracle. The
  mid-2020s mmCIF-as-default transition broke many parsers in subtle
  ways; pinning two independent oracles caught at least one such
  regression during proteon's pdbtbx upgrade work.
- Tolerances are a load-bearing part of the claim, not an
  engineering convenience. The 0.001 Å bound is the PDB format
  precision floor — anything stricter would fail by physics, anything
  looser would let real coord drift hide.
