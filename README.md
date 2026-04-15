# Ferritin

**A structural bioinformatics compute kernel** — the hot inner loop for loading, analyzing, and aligning protein structures. Rust at the core, Python and CLI on top. You bring your data, your pipeline, your infra; ferritin brings the compute.

- **Compute kernel, not a platform.** Outputs are numpy, Arrow, and Parquet. No daemons, no services, no opinions about your storage or your scheduler. You call it; it returns data.
- **Rust core, multiple frontends.** One kernel, exposed as a Rust library, a Python package, CLI binaries, and an Arrow/Parquet ingester. Python is one surface, not the product.
- **Validated under load.** 45 100 real PDBs processed end-to-end in under 18 minutes on 120 cores; oracle-tested against BALL, USAlign, and Biopython; 10 correctness bugs found and fixed by running at scale.

```python
import ferritin

s = ferritin.load("1crn.pdb")

# Analyze
sasa        = ferritin.atom_sasa(s)              # Shrake-Rupley
ss          = ferritin.dssp(s)                   # Kabsch-Sander
phi, psi, _ = ferritin.backbone_dihedrals(s)     # Ramachandran
hbonds      = ferritin.backbone_hbonds(s)        # H-bond detection

# Align
r = ferritin.tm_align(s1, s2)
print(f"TM-score: {r.tm_score_chain1:.4f}")

# Batch: thousands of structures, all cores, zero GIL
results = ferritin.load_and_analyze(pdb_files, n_threads=-1)
```

## Installation

```bash
pip install ferritin
```

## Validated Under Load

Ferritin is battle-tested against real data and independent oracles, not just happy-path unit tests.

### 50K real PDB corpus (120 cores, 250 GB RAM)

| Stage | Throughput | Notes |
|-------|-----------|-------|
| Load (PDB + mmCIF) | 45 100 / 50 000 in ~7 min | 4 900 filtered by size / atom count |
| SASA | 160 structures/s | Shrake-Rupley, ProtOr radii |
| DSSP | 2 810 structures/s | Kabsch-Sander |
| Backbone dihedrals | 3 612 structures/s | — |
| Backbone H-bonds | 3 180 structures/s | — |
| Hydrogen placement | 335 structures/s | 25.2 M hydrogens placed |
| Prepare (H + LBFGS minimize) | 4.8 structures/s | 152 / 200 converged in 200 LBFGS steps |
| **End-to-end** | **45 100 structures in 17.9 min** | full pipeline, single run |

Running at 50K caught 10 correctness bugs that no smaller-scale test had surfaced — OOM paths on NMR ensembles, analytical force sign errors, neighbor-list grid overflows, switch-function derivative mistakes. Each was fixed with a regression test.

### Independent oracles

| Check | Baseline | Agreement |
|-------|----------|-----------|
| TM-align | C++ USAlign | 0.003 median TM-score drift over 4 656 pairs |
| SASA | Biopython | within 0.17 % of Biopython median on 1 000 PDBs |
| AMBER96 energy | BALL (C++ force field library) | 8 / 8 crambin components match within 0.1 % |
| PDB parsing | Random 1 000 PDB archive subset | 96.3 % load success |

More oracles planned (FreeSASA, DSSP-3, PDBFixer, OpenMM) on a curated reference set.

## What's in the kernel

| Capability | Entry point |
|------------|-------------|
| **Structure I/O** (PDB, mmCIF) | `ferritin.load("protein.pdb")`, `ferritin.save(s, "out.pdb")` |
| **TM-align / US-align family** | `tm_align`, `soi_align`, `flex_align`, `mm_align` (multi-chain) |
| **SASA** (Shrake-Rupley, ProtOr) | `atom_sasa`, `residue_sasa`, `relative_sasa`, `total_sasa` |
| **DSSP** (Kabsch-Sander) | `dssp(s)` — no external binary |
| **Backbone dihedrals** | `phi, psi, omega = backbone_dihedrals(s)` |
| **H-bond detection** | `backbone_hbonds(s)`, `geometric_hbonds(s)` |
| **Force field + minimize** | AMBER96 / CHARMM19-EEF1, LBFGS minimizer, `prepare()` pipeline |
| **Atom selection language** | `select(s, "backbone and resid 1-50 and chain A")` |
| **Geometry & analysis** | `contact_map`, `distance_matrix`, `radius_of_gyration`, `kabsch_superpose` |
| **DataFrame export** | `to_dataframe(s, engine="pandas" \| "polars")` |

Every capability has a batch variant (`batch_dssp`, `batch_total_sasa`, `batch_prepare`, ...) that runs across structures on all cores with the GIL released.

Short runnable scripts for each area live in [`examples/`](examples/).

## Integration Surface

Ferritin is a library, not a service. You pick how to feed it data and where the outputs go.

**Python** — `pip install ferritin`, import, call. All Python functions release the GIL on the hot path.

**CLI binaries** — structure alignment and bulk ingestion:

```bash
ferritin-align s1.pdb s2.pdb                             # TM-align / US-align
ferritin-align chain_list.txt --dir /pdbs/ --outfmt 2    # all-vs-all, parallel
ferritin-ingest structures/ --out features.parquet       # bulk → Parquet
```

**Arrow / Parquet** — zero-copy columnar output, ready for DuckDB, polars, pandas, Spark, or PyTorch Geometric:

```python
import ferritin

s = ferritin.load("protein.pdb")

ipc = ferritin.to_arrow(s, "1crn")               # Arrow IPC bytes (zero copy)
ferritin.to_parquet(s, "out.parquet", "1crn")    # direct-to-disk, Zstd
structs = ferritin.from_parquet("out.parquet")   # round-trip
```

```python
# Query a bulk-ingested Parquet directly — no loader code needed
import duckdb
duckdb.sql("""
    SELECT structure_id, residue_serial, b_factor
    FROM 'features.parquet'
    WHERE residue_name = 'GLY' AND atom_name = 'CA'
      AND chain_id = 'A' AND b_factor > 30
    ORDER BY b_factor DESC
""")
```

**Rust** — depend on the crates directly. See Architecture.

## Architecture

```
Pure Rust (no Python dependency)
├── ferritin-align     TM-align, US-align, SOI-align, FlexAlign, MM-align
│   ├── core/          Kabsch, TM-score, DP, secondary structure
│   └── ext/           US-align extensions (BLOSUM, SOI, Flex, MM-align)
├── ferritin-io        PDB / mmCIF I/O via pdbtbx
└── ferritin-bin       CLI binaries

PyO3 connector (cdylib, rayon, GIL-released)
└── ferritin-connector SASA, DSSP, H-bonds, force field, alignment wrappers

Python package
└── ferritin           Pythonic API over the connector
```

## Agent-Aware Documentation

Ferritin includes structured **Agent Notes** on selected public-boundary functions — short guidance for AI agents and LLM tools that call ferritin programmatically.

```python
def atom_sasa(structure, probe=1.4, n_points=960):
    """Compute per-atom SASA using Shrake-Rupley.

    Agent Notes:
        WATCH: probe=1.4 is the standard water probe radius. Changing it
               changes the physical interpretation of the result.
        PREFER: For many structures, use batch_total_sasa() with n_threads=-1.
        COST: Crambin (327 atoms): ~12 ms. Large complex (58 k atoms): ~230 ms.
    """
```

| Prefix | Meaning |
|--------|---------|
| `WATCH` | Common gotchas that are not bugs |
| `PREFER` | Better alternatives for common use cases |
| `COST` | Performance or memory implications |
| `INVARIANT` | Guarantees callers can rely on |

Policy: used only on public APIs where misuse is likely or where a caller must make a decision. Max three notes per function, no duplication with the docstring.

## Examples

See [`examples/`](examples/) for runnable scripts — I/O, alignment, contact maps, Ramachandran analysis, SASA, and geometric-deep-learning pipelines (the latter use ferritin for data prep plus PyTorch Geometric for training; ferritin itself has no DL dependency).

## References

**Alignment**
- Zhang & Skolnick. "TM-align: a protein structure alignment algorithm based on the TM-score." [*Nucleic Acids Research* 33, 2302-9 (2005)](https://doi.org/10.1093/nar/gki524)
- Zhang, Pagnon Braunstein, et al. "US-align: universal structure alignments of proteins, nucleic acids, and macromolecular complexes." [*Nature Methods* 19(9), 1109-1115 (2022)](https://doi.org/10.1038/s41592-022-01585-1)

**Structural analysis**
- Kabsch & Sander. "Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features." [*Biopolymers* 22, 2577-2637 (1983)](https://doi.org/10.1002/bip.360221211)
- Shrake & Rupley. "Environment and exposure to solvent of protein atoms." [*J Mol Biol* 79(2), 351-71 (1973)](https://doi.org/10.1016/0022-2836(73)90011-9)
- Tien et al. "Maximum allowed solvent accessibilities of residues in proteins." [*PLoS ONE* 8(11), e80635 (2013)](https://doi.org/10.1371/journal.pone.0080635)

**Force field**
- Cornell et al. "A Second Generation Force Field for the Simulation of Proteins, Nucleic Acids, and Organic Molecules." [*J. Am. Chem. Soc.* 117(19), 5179-5197 (1995)](https://doi.org/10.1021/ja00124a002)
- Lazaridis & Karplus. "Effective energy function for proteins in solution." [*Proteins* 35(2), 133-152 (1999)](https://doi.org/10.1002/(SICI)1097-0134(19990501)35:2%3C133::AID-PROT1%3E3.0.CO;2-N) — CHARMM19-EEF1 implicit solvent

**Infrastructure**
- Hildebrandt et al. "BALL — Biochemical Algorithms Library 1.3." [*BMC Bioinformatics* 11, 531 (2010)](https://doi.org/10.1186/1471-2105-11-531)
- Schulte, D. "pdbtbx: A Rust library for reading, editing, and saving crystallographic PDB/mmCIF files." [*JOSS* 7(77), 4377 (2022)](https://doi.org/10.21105/joss.04377)

## License

MIT — see [LICENSE](LICENSE). Clean-room Rust reimplementation. Please cite the papers referenced above if you use ferritin in published work.
