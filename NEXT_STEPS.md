# NEXT_STEPS.md

Pick-up notes after the v0.1.1 publish session (2026-04-18). Ordered by
urgency, not importance — "do first" is at the top.

---

## Right now — finish the v0.1.1 release

**The macos-13 x86_64 wheel job has been queued for ~2+ hours** and is
blocking the PyPI publish step. Either decision is fine:

- **Cancel it.** Intel Mac is dropped from the matrix for v0.1.2+ anyway
  (see `.github/workflows/release.yml` comment). v0.1.1 ships with 4 wheel
  targets + sdist, Intel Mac users wait one version.
  ```bash
  gh run cancel $(gh run list --workflow=release.yml --limit 1 --json databaseId --jq '.[0].databaseId')
  ```
  The publish jobs will then resume using `fail-fast: false` — watch:
  ```bash
  gh run watch $(gh run list --workflow=release.yml --limit 1 --json databaseId --jq '.[0].databaseId')
  ```
- **Wait it out.** macos-13 runner capacity is unpredictable — could land
  in 10 min, could sit another day.

**After publish succeeds**, verify end-to-end in a fresh venv:
```bash
python3.12 -m venv /tmp/ferritin-smoke && source /tmp/ferritin-smoke/bin/activate
pip install ferritin
python -c "import ferritin; s = ferritin.load('test-pdbs/1crn.pdb'); print(s.atom_count)"
```
A green CI only proves wheels built, not that they import.

---

## This week

### 1. Release announcement

First PyPI release is a moment; a package released to zero announcement
stays at zero users. Short post for HN / r/bioinformatics / Mastodon / lab
Slack — what ferritin does, validation story (0.2% AMBER96 parity, 99.1%
50K-PDB battle test, 0.003 TM-score drift), 10-line code example.
Memory: this was the user-owned #1 priority before the citation detour.

### 2. Upstream outreach (before the announcement)

Short heads-up emails to the authors whose work ferritin ports, **before**
anyone else hears about it on HN. One paragraph each, linking the repo
and the specific validation claim that lands on their tool. In rough
priority:

- **Martin Steinegger** (MMseqs2, Foldseek) — most load-bearing; our
  search layer is a port, the GPU kernel design follows libmarv, and
  he's an author on both citations.
- **Yang Zhang lab** (TM-align, US-align) — alignment core.
- **Felix Kallenborn** (libmarv / GPU-MMseqs2) — we cite the paper and
  name the kernel "libmarv".
- **Peter Eastman** (OpenMM) — our primary force-field oracle.
- **Andreas Hildebrandt** (BALL) — CHARMM19 oracle.

See `tests/oracle/README.md` for install pointers and `docs/WHY.md §"On
Credit and Invisible Debt"` for the framing.

### 3. Zenodo ↔ GitHub integration

User-side UI flip:
1. Log in to <https://zenodo.org> (GitHub auth).
2. Settings → GitHub → flip the switch for `theGreatHerrLebert/ferritin`.
3. The next GitHub release (v0.1.2) auto-archives with a DOI.
4. Add the DOI to `CITATION.cff` as `identifiers`.

Purely user action, no repo change needed until the DOI shows up.

### 4. Cut v0.1.2 cleanly

When v0.1.1 is settled, next release takes advantage of the matrix
simplifications this session landed:

- Bump `Cargo.toml` `[workspace.package].version` + both `pyproject.toml`
  versions to `0.1.1` → `0.1.2`.
- Bump `packages/ferritin/pyproject.toml` runtime dep to
  `ferritin-connector>=0.1.2`.
- Tag `v0.1.2`, create release, workflow fires — no macos-13 job, so the
  whole run should land in < 20 min.

---

## Next session — small engineering tasks

### Reduce oracle (the one we signaled intent on)

Design questions already captured in `tests/oracle/README.md` §Intent:

1. Ferritin polar-only H for CHARMM19 vs reduce's full H set →
   parametrize over FF.
2. Asn/Gln/His flip search → `-NOFLIP` or filter flip-dependent atoms.
3. Tolerance: per-H coordinate RMSD, default to that.

Once those are resolved, ~100-line `test_reduce_hydrogen_oracle.py` +
promote from candidates table to install table.

### Second-wave candidate oracles

Listed in `tests/oracle/README.md` §Candidate. Order of easiest-first:

- **DSSP binary (`mkdssp`)** — probably a few hours. Ferritin ships its
  own DSSP port; per-residue H/E/G/I agreement with upstream `mkdssp` is
  a trivial oracle. Low risk, high coverage.
- **PDB2PQR** — protonation / pKa / charges. Medium effort, only useful
  when users hit non-default-pH electrostatics.
- **MolProbity** — full validation suite (bundles reduce + probe). Heavy
  install, pays off once ferritin-prepped structures are used at scale.

### Slow test marker

We added `@pytest.mark.oracle("tool")` this session; complement with
`@pytest.mark.slow` for tests >10s, so dev loops can
`pytest -m "not slow and not oracle"` and get a fast feedback cycle.
Register in `tests/conftest.py` alongside the oracle marker.

---

## Longer arc

### crates.io publish

Blocked: `pdbtbx` is a git dep (`theGreatHerrLebert/pdbtbx @ c82e8c0`),
and crates.io rejects git deps. Two paths out:

1. Upstream the needed pdbtbx patches to `douweschulte/pdbtbx`, wait for
   a crates.io release that includes them, then `cargo publish` the
   ferritin crates. Real engineering commitment (multi-week, depends on
   upstream responsiveness).
2. Fork-publish as `pdbtbx-ferritin` on crates.io. Permanent maintenance
   tax.

Neither is urgent — PyPI is the distribution channel that matters for the
primary (Python) audience. Revisit when there's explicit crates.io demand.

### Foldseek alphabet — close the recall gap

Current state: ~15% recall gap vs upstream Foldseek at TM≥0.5,
near-parity at TM≥0.9. Benchmarks at
`validation/foldseek_ferritin_5k_50q_union50.report.md`. The training
scripts are `validation/train_vqvae.py` and `validation/train_alphabet.py`.
If/when this closes, Foldseek moves from "experimental alphabet inspired
by" to "faithful port with X% parity" in the README.

### JOSS paper

When v0.2+ has feature stability (more oracle tests wired, reduce + DSSP
+ PDB2PQR in, the alphabet gap closed or explicitly scoped), the project
is at JOSS-submission maturity. Paper body largely mirrors `docs/WHY.md`
+ the validation numbers. Multi-week effort, but gives ferritin its own
citable reference instead of piggybacking on the tools it ports.

### SOTA science + GPU optimizations

Per the last session handoff (`devdocs/NEXT_SESSION.md`). Not blocked by
anything this session produced; resume whenever.

---

## Session audit trail (what this session produced)

For continuity when picking this up from memory alone:

- v0.1.0 manually uploaded to PyPI to claim the project names.
- Trusted publishers configured for both `ferritin` and
  `ferritin-connector` (latter via project-settings page after the pending
  UI crashed).
- v0.1.1 tagged + release created; workflow fixed twice (maturin
  interpreter flag; macos-13 pyo3 build cap).
- Dropped Python 3.10/3.11 (matrix was 8 jobs, now 4).
- Dropped Intel Mac from Rust test matrix (permanently) and from release
  matrix (v0.1.2+).
- Oracle pattern made first-class: `@pytest.mark.oracle` marker,
  `tests/oracle/README.md` with install pointers + candidate table,
  `CONTRIBUTING.md`, `CITATION.cff`, Acknowledgements paragraph, WHY.md
  "On Credit and Invisible Debt" section.
- Per-file citations added to 15 algorithm modules (Kabsch, TM-align,
  US-align, MMseqs2 stages, BLOSUM/PAM/VTML, CHARMM19/EEF1/AMBER96/OBC GB,
  SASA, DSSP, hbond, libmarv GPU SW).
- Foldseek + libmarv papers cited, module-level and README-level.
