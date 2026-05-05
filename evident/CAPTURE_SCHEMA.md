# EVIDENT capture schema

This document defines what each oracle runner captures per PDB when it
runs, the on-disk layout the captures land in, and how downstream
renderers + documentation consume them. The intent: every release-tier
sweep produces a substrate rich enough to derive HTML reports,
paper-ready tables, outlier walk-throughs, and cross-release
comparisons **without re-running the compute**.

The schema is versioned. Current version: **v1**. Bumps require a
deprecation window where renderers must accept both old and new shapes
for at least one release.

## TL;DR

- Three capture tiers: `minimal`, `extended`, `full`. Selected via env
  var `PROTEON_CAPTURE_LEVEL`.
- One **JSONL summary file** per claim per run. One record per PDB.
- Adjacent **`hardware.json`** captures host fingerprint once per run.
- For `full`-tier runs (or `extended`-tier runs that opt-in), per-PDB
  artifacts (intermediate structures, trajectories, per-residue arrays)
  land in a sibling `pdbs/<pdb_id>/` directory.
- `lock_release_replays.py` sha256-pins the summary JSONL +
  `hardware.json`. The per-PDB tree is referenced by relative path; its
  hash chain is opt-in (rolled into the summary JSONL via per-record
  `pdbs_dir_sha256` if `--with-pdbs-tree` is passed at lock time).

## Three tiers

```
minimal   : headline floats + status, ~500 B/PDB
extended  : + per-residue / per-stage / hardware fingerprint, ~5 KB/PDB
full      : + minimization trajectories + intermediate PDBs, ~50 KB/PDB
```

| Tier | When to use | 50K total disk |
|---|---|---|
| `minimal` | CI gates, fast regression guards. Default. | ~25 MB |
| `extended` | Release-tier sweeps. Default for v0.2.0 work. | ~250 MB |
| `full` | Forensic analysis, paper figures, outlier deep-dives, smaller corpora (≤1k). | ~2.5 GB at 50K |

Tier is opt-in via env var:

```bash
PROTEON_CAPTURE_LEVEL=extended python validation/charmm19_eef1_ball_oracle.py
```

Default is `minimal` so existing runners stay backward-compatible.
Setting `PROTEON_CAPTURE_LEVEL=full` implicitly enables `extended`.

## Filesystem layout

For a release sweep, each claim's runner writes into a single
directory:

```
validation/v0.2.0_runs/
  charmm19_corpus_50k/                            # one dir per claim
    summary.jsonl                                  # 1 record per PDB
    hardware.json                                  # host fingerprint, written once
    schema_version: 1
    pdbs/                                          # only at extended/full tier
      1abc.pdb/                                    # per-PDB drill-down dir
        input_polh.pdb                             # the input we evaluated
        proteon_components.json                    # per-component proteon energies
        ball_components.json                       # per-component oracle energies (if applicable)
        per_residue.npz                            # extended only, if claim has per-residue data
        minimization_trace.npz                     # full only, if claim minimizes
        meta.json                                  # per-PDB hardware + wall-time breakdown
```

The directory name `<release-tag>_runs/<claim-key>/` is the
**capture directory**. Its layout is the contract; the runner picks the
exact path. `claim-key` is a short slug (e.g. `charmm19_corpus_50k`).

`summary.jsonl` is the canonical artifact — `lock_release_replays.py`
sha256-pins it as the claim's evidence file. The `pdbs/` tree is
**referenced** by summary records (each `summary.jsonl` line includes
a `pdb_artifact_dir: "pdbs/1abc.pdb"` field when present), but isn't
sha256-pinned by default unless `--with-pdbs-tree` is passed.

## Per-record JSONL schema

### `minimal` tier (current behaviour, preserved)

```jsonc
{
  "schema_version": 1,
  "schema_tier": "minimal",
  "pdb": "1abc.pdb",
  "status": "ok" | "skipped" | "error" | "load_error" | "pipeline_error",
  "skipped": false,                              // optional, true when status == skipped
  "error": "...",                                // present iff status indicates failure
  "wall_s": 0.42                                 // total wall time
}
```

Plus claim-specific headline floats. Examples:

- corpus oracle: `proteon`, `ball`, `rel_diff` dicts (per-component)
- 50K battle test: `initial_energy`, `final_energy`, `steps`,
  `converged`, `atoms`
- SASA: `biopython_total`, `proteon_total`, `relative_diff`,
  `n_atoms`, `speedup`

### `extended` tier (adds to `minimal`)

```jsonc
{
  // ... all minimal fields, plus:
  "schema_tier": "extended",
  "structure_meta": {
    "n_atoms": 1413,
    "n_residues": 154,
    "n_chains": 1,
    "n_disulfides": 0,
    "n_HIS": 4,                                  // helps interpret CHARMM outliers
    "ss_helix_frac": 0.71,                       // from DSSP, helps interpret SASA outliers
    "ss_strand_frac": 0.05,
    "max_residue_idx_gap": 0                     // detects PDB chain breaks
  },
  "wall_breakdown_s": {
    "load": 0.05,
    "prep_polar_h": 0.12,
    "proteon_compute": 0.08,
    "oracle_compute": 0.31,                      // ball/openmm/biopython side
    "io": 0.02
  },
  "diagnostics": {
    "gradient_norm_kJ_per_mol_per_A": 142.7,    // non-zero on un-relaxed inputs
    "rmsd_input_minimized_A": 0.31,             // post-minimization RMSD vs input, if applicable
    "n_unassigned_atoms": 0
  },
  "hardware_record_id": "monster3-rtx5090-w0",  // joins to hardware.json
  "pdb_artifact_dir": "pdbs/1abc.pdb"           // null if no drill-down captured
}
```

### `full` tier (adds to `extended`)

The summary JSONL gains pointers; the actual data lives under
`pdbs/<pdb>/`:

```jsonc
{
  // ... all extended fields, plus:
  "schema_tier": "full",
  "trace": {                                     // null if claim doesn't minimize/MD
    "captured_steps": 50,
    "trace_file": "pdbs/1abc.pdb/minimization_trace.npz"
  },
  "intermediate_structures": [
    "pdbs/1abc.pdb/input_polh.pdb",
    "pdbs/1abc.pdb/minimized.pdb"
  ]
}
```

The `.npz` for trajectories has shape:

```python
np.load("minimization_trace.npz")
# → {
#     "step":          int32[N],            # step index
#     "energy_kJmol":  float64[N],          # total energy
#     "grad_norm":     float64[N],          # gradient norm
#     "rmsd_to_init":  float64[N],          # if available, else absent
#   }
```

Per-residue arrays for SASA and similar:

```python
np.load("per_residue.npz")
# → {
#     "residue_id":      str[R],            # "A:42:LYS"-style
#     "proteon_value":   float64[R],
#     "oracle_value":    float64[R],
#     "rel_diff":        float64[R]
#   }
```

## `hardware.json`

Written once per claim per run. Records the host environment in detail
beyond what the lock-time `pip freeze` captures:

```jsonc
{
  "schema_version": 1,
  "captured_at_utc": "2026-05-04T08:30:00Z",
  "host": "monster3",
  "os": "AlmaLinux 10.1 (Linux 5.15.0-...)",
  "cpu": {
    "model": "AMD Threadripper PRO 7975WX",
    "logical_cores": 64,
    "physical_cores": 32
  },
  "ram_gb": 256,
  "gpu": {
    "model": "NVIDIA GeForce RTX 5090",
    "vram_mb": 32607,
    "driver": "595.58.03",
    "cuda_runtime": "12.8",
    "compute_capability": "9.0"
  },
  "storage_backend": "/globalscratch on hoard1 (NFS)",
  "container_runtime": "apptainer 1.4.0",
  "image_digest": "sha256:abc123...",          // null if running native venv
  "image_tag": "ghcr.io/...:cuda-v0.1.4",       // null if running native venv
  "n_workers": 32,
  "task_timeout_s": 120,
  "captured_records": [                         // hardware_record_id values used in summary.jsonl
    "monster3-rtx5090-w0"
  ]
}
```

`hardware_record_id` exists per record because some sweeps may span
multiple hosts or GPUs (e.g. a multi-host fold-preservation run); the
join key keeps the per-record audit clean. For single-host runs there
is exactly one entry.

## Renderer contract

Each per-claim renderer (`validation/report/render_*.py`) declares
what tier it requires. Rendering against a record from a lower tier
should degrade gracefully (skip the section, do not crash):

```python
# validation/report/render_corpus_oracle.py
TIER_REQUIRED = "minimal"
TIER_PREFERRED = "extended"

def render(claim, summary_jsonl, hardware_json=None):
    records = load_jsonl(summary_jsonl)
    if records[0].get("schema_tier") == "extended":
        # render the per-stage wall breakdown plot
        ...
    if hardware_json:
        # render the hardware fingerprint section in the header
        ...
```

## Lock-time integration

`evident/scripts/lock_release_replays.py`:

1. Reads the claim's `evidence.artifact` path. For schema-v1 claims
   pointing at a directory rather than a single file, treats the
   directory's `summary.jsonl` as the canonical artifact and the
   directory's `hardware.json` as a sibling.
2. sha256-pins **both** files into the manifest's claim record.
3. If the directory has a `pdbs/` subdirectory, records its existence
   + a top-level sha256 of the directory's structure
   (`find pdbs/ -type f | sort | xargs sha256sum | sha256sum`) but
   does NOT include individual file hashes by default. Pass
   `--with-pdbs-tree` to enable per-file hashing (slow on full-tier
   sweeps; ~50K files).
4. The renderer is invoked with both files; it picks the level it
   can support.

Backward compatibility: claims pointing at a single `.jsonl` file
(today's format) continue to work — the locker treats them as
schema-v0 / `minimal`-tier and the existing renderers consume them
unchanged.

## Documentation derivation pipeline

Given a release with `extended`-tier capture, four kinds of
documentation can be auto-generated:

### 1. Per-claim HTML reports (already exist for 4 claims)

Renderers gain new sections when records carry `extended` data:

- **Per-stage wall-time breakdown** — diagnoses where time goes
- **Hardware fingerprint** — banner on each per-claim HTML
- **Outlier metadata table** — pulls `structure_meta` for the top-N
  worst-rel-diff PDBs, surfacing patterns ("4 of the 5 worst CHARMM
  outliers have ≥3 disulfides")
- **Per-residue heatmaps** — for SASA/DSSP claims with per-residue
  arrays in `pdbs/<pdb>/per_residue.npz`

### 2. Outlier walk-through markdown

`evident/scripts/render_outliers.py <release-tag> <claim-id> --top-n=10`
emits one markdown file per outlier PDB with:

- structural metadata
- the rel-diff numbers
- the link into the per-PDB artifact dir
- a one-line interpretive comment if a matching pattern is registered
  (e.g. "free CYS HG residue-template miss")

The markdowns land at `evident/reports/<tag>/outliers/<claim>/<pdb>.md`
and link from the per-claim HTML.

### 3. Paper / poster ready tables

`evident/scripts/export_tables.py <release-tag> --format=csv|tex`
emits per-claim CSV/TeX tables with the headline metrics. Schema:

```
pdb,n_atoms,proteon_total,oracle_total,rel_diff,wall_s
1abc.pdb,1413,-12345.6,-12350.2,0.000037,0.42
...
```

### 4. Cross-release comparison plots

`evident/scripts/compare_releases.py v0.1.3 v0.1.4 --claim=charmm19-corpus-50k`
overlays the relative-diff distributions from two release bundles. The
output PNG goes into `evident/reports/<latest-tag>/comparisons/`.

## Storage and publishing

- **summary.jsonl** ≤ 250 MB per release at extended tier ⇒
  comfortable to commit to the repo (gitignored on the source tree;
  bundled into release artifacts and Pages).
- **hardware.json** is small (~2 KB per run); always committed.
- **pdbs/** tree at full tier reaches multi-GB; **publish via GH
  Release assets**, not the repo. The lock manifest references the
  asset URL; sha256 is recorded for the directory bundle as a tarball.
- **Pages** serves the per-claim HTML from
  `evident/reports/<tag>/`; paths into `pdbs/<pdb>/` resolve to
  external URLs (release asset download links) when the user clicks
  through.

## Running in containers (v0.2.0 data-mount contract)

The published EVIDENT image (`ghcr.io/thegreatherrlebert/proteon-evident:<tag>`)
bundles every reference oracle (`ball-py`, OpenMM, mkdssp, USAlign, GROMACS,
MMseqs2, …) and every release-tier runner script. **The image deliberately
does NOT bundle the input corpus** — corpus tarballs balloon image size,
expire faster than the image (PDB sample SHAs change as wwPDB grows), and
make the same image less reusable across replay scenarios.

Instead, the image exposes a **bind-mount contract**:

| Mount | Role | Env var the entrypoint sets |
|---|---|---|
| `/data/pdbs` | Read-only input corpus directory (PDBs/CIFs) | `PROTEON_CORPUS_DIR` |
| `/data/out`  | Read-write output directory (per-claim JSONLs)  | `PROTEON_OUTPUT_DIR` |

The image entrypoint (`evident/scripts/docker-entrypoint.sh`) checks for
these well-known mountpoints at startup. If present, it exports the
matching env vars before dispatching. If absent, it unsets them — runners
fall back to their repo-relative defaults (so a smoke run with no mounts
still does something deterministic on the bundled `tests/fixtures/` data).

User-supplied `-e PROTEON_CORPUS_DIR=...` overrides the convention:

```bash
# Convention (recommended):
docker run --rm \
  -v $(pwd)/my_pdbs:/data/pdbs \
  -v $(pwd)/my_out:/data/out \
  ghcr.io/thegreatherrlebert/proteon-evident:v0.2.0 \
  replay proteon-charmm19-fold-preservation-vs-openmm-release-1k-pdbs

# Explicit (when paths can't follow /data/* — e.g. SLURM / Apptainer):
apptainer exec \
  --bind /scratch/dteschner/pdbs:/corpus \
  --bind /scratch/dteschner/out:/results \
  --env PROTEON_CORPUS_DIR=/corpus \
  --env PROTEON_OUTPUT_DIR=/results \
  proteon-evident_v0.2.0.sif \
  replay <claim-id>
```

### Runner contract

Every release-tier runner under `validation/` honors this precedence:

1. Legacy per-runner env var if defined (e.g. `PROTEON_PDB_DIR`,
   `PROTEON_CHARMM_ORACLE_OUT`) — keeps existing monster3 batch scripts
   working unchanged.
2. v0.2.0 universal env vars `PROTEON_CORPUS_DIR` / `PROTEON_OUTPUT_DIR`.
3. Repo-relative fallback (`validation/pdbs_1k_sample/`,
   `validation/<runner-canonical-filename>.jsonl`).

Each runner picks a canonical output filename within `PROTEON_OUTPUT_DIR`,
documented in its module docstring. The joiner
(`validation/fold_preservation/join_fold_preservation.py`) reads the same
directory the runners wrote to (also via `PROTEON_OUTPUT_DIR`), so the
chain stays consistent.

### Auxiliary tool paths

For oracles that ship as vendored binaries inside the image, the
entrypoint also exports a small set of well-known tool paths so runners
can locate them without hardcoding `/usr/local/bin/...`:

| Env var | Image default | Used by |
|---|---|---|
| `USALIGN_BIN` | `/usr/local/bin/USalign` | `validation/run_validation.py` (SASA-vs-Biopython) |

User overrides take precedence in the same way as the data-mount vars.

### What this contract is NOT

- **Not a configuration system.** Runners take exactly two inputs from
  the environment: a corpus directory and an output directory. Anything
  else (worker count, seed, claim-specific knobs) stays in the runner's
  own argparse / env vars.
- **Not a substitute for the source-tree replay.** A reviewer who wants
  to step through the Python source can still `git checkout v0.2.0 &&
  pip install -e . && python validation/...py` exactly as before. The
  image is for "I want a verified result, not a development environment".

## Migration path

Per-claim adoption rolls out incrementally:

1. **Schema doc landed (this PR)** — formal spec exists, no runner
   changes yet.
2. **Reference implementation: corpus oracle** —
   `validation/charmm19_eef1_ball_oracle.py` gains
   `PROTEON_CAPTURE_LEVEL` env var. Default `minimal` (no behaviour
   change). Other tiers add the new fields.
3. **Renderer extension** — `render_corpus_oracle.py` learns to read
   `extended`-tier records. Renders the new sections when present;
   silently skips when records are minimal.
4. **Other runners follow** — `amber96_oracle.py`, `run_validation.py`
   (SASA), `tm_fold_preservation_*.py`, etc. Each PR is one runner +
   any required renderer changes.
5. **`lock_release_replays.py` directory-mode** — once at least one
   claim has `extended`-tier artifacts, teach the locker to handle
   the directory + `hardware.json` + optional `pdbs/` layout.

A claim is "extended-tier ready" when its runner emits the schema, its
renderer reads the schema, and a sample run on a small corpus
validates end-to-end. Once all release-tier claims are
extended-tier ready, **v0.2.0 captures everything in one sweep**.

## Open questions / future work

These are not blockers for v0.2.0 capture but should be revisited as
the schema matures:

- **MD trajectory capture**: when proteon ships a real MD claim
  (Langevin / Verlet, not just minimization), the schema needs a
  `traj` block with frame stride, energy / temperature / RMSD
  trajectories at the chosen stride. Spec'd loosely above; concrete
  fields TBD when the first MD claim lands.
- **Force-field component decomposition for minimized structures**:
  currently we capture initial + final energy, not per-component.
  Worth adding for `extended` tier since it's cheap and unlocks
  "which component is dominating the gradient" questions.
- **Cross-claim joins**: a structure that fails the corpus oracle
  AND fails fold-preservation likely has a common root cause. A
  small SQL-style query layer (`evident/scripts/cross_claim_query.py`)
  would let us answer "which PDBs fail in both X and Y?" without
  hand-joining JSONLs.
- **Provenance of capture artifacts**: today the runner writes
  whatever it computes. Future tightening: each captured PDB-dir
  carries a small per-pdb manifest (`meta.json`) with sha256 of
  every file in that dir, so a reviewer can verify a single
  drill-down without trusting the parent summary.
