# Ferritin Roadmap

**Last updated: 2026-04-08**

Ferritin is a pipeline-native, high-performance compute layer for structural bioinformatics. Rust core, Python bindings, batch-first design.

---

## Done

### Core Pipeline (production-grade, 5K validated)
- [x] I/O: PDB/mmCIF load/save, batch parallel, tolerant loading
- [x] SASA: Shrake-Rupley, Bondi + ProtOr radii (0.01% match to FreeSASA)
- [x] DSSP: Kabsch-Sander, uses real H when present
- [x] H-bonds: backbone (Kabsch-Sander energy) + geometric (distance-based)
- [x] Dihedrals: phi/psi/omega with backbone-break detection
- [x] Selection language: MDAnalysis/PyMOL-style boolean queries
- [x] Alignment: TM-align, US-align, SOI-align, FlexAlign, MM-align
- [x] Arrow/Parquet export: zero-copy interop with pandas/polars/DuckDB
- [x] Batch parallelism: every feature has `batch_*` + `load_and_*` variants

### Hydrogen Placement (BALL-validated)
- [x] Phase 1: backbone amide N-H (DSSP bisector, 1.01 Å)
- [x] Phase 2: sidechain templates for all 20 standard amino acids
- [x] Phase 3: general algorithm (bond orders, ring detection, MMFF94 bond lengths)
- [x] Disulfide detection: skip CYS HG when SG-SG < 2.5 Å
- [x] Terminal H: N-terminal NH3+ (H1/H2/H3), C-terminal HXT
- [x] Water H: optional `include_water=True`
- [x] Oracle-validated: exact match with BALL (316=316 on crambin)

### Missing Atom Reconstruction
- [x] Fragment templates for 20 standard amino acids (from BALL FragmentDB)
- [x] 3-point rigid body superposition (match_points)
- [x] BFS over template bond graph for coordinate propagation
- [x] reconstruct_fragments() adds missing heavy atoms + H

### Force Field & Simulation
- [x] AMBER96 energy: bond stretch, angle bend, torsion, LJ 12-6, Coulomb
- [x] Steepest descent minimizer with adaptive step size
- [x] Conjugate gradient minimizer (Polak-Ribiere + Armijo line search)
- [x] Velocity Verlet MD integrator (NVE)
- [x] Berendsen thermostat (NVT)
- [x] Bond order estimation from distances
- [x] Ring atom detection (BFS cycle finding)

### Infrastructure
- [x] Parallel validation runner (batch-first, FreeSASA + Gemmi + PDBFixer oracles)
- [x] 5K structure validation (100% load, 99.6% DSSP, 100% H-bonds)
- [x] Agent notes: condensed PREFER/COST/WATCH style

---

## Next Up

### Validation & Testing (highest priority)
- [ ] Re-run 5K validation with bond length fix + disulfide fix + ProtOr radii
- [ ] BALL oracle comparison on 50+ structures for Phase 3 general H placement
- [ ] DSSP oracle: wire GROMACS `gmx dssp` or install mkdssp
- [ ] CG minimizer convergence comparison vs steepest descent
- [ ] Bond order accuracy validation against known structures
- [ ] Missing atom reconstruction validation on structures with incomplete residues
- [ ] MD energy conservation test (longer NVE runs, drift analysis)
- [ ] Expand Gemmi monomer library for full CCD coverage

### Structure Preparation Gaps (vs PDBFixer)
- [ ] Protonation state assignment (pH-dependent His tautomers, Asp/Glu)
- [ ] Asn/Gln/His flip optimization (amide/ring orientation)
- [ ] Fragment variant selection (N-terminal, C-terminal, disulfide CYS variants)
- [ ] Nucleotide support in fragment templates
- [ ] Non-standard residue templates (MSE, PTR, HYP, etc.)

### Force Field Improvements
- [ ] Torsion gradient (currently energy-only, no forces for torsions)
- [ ] L-BFGS minimizer (faster convergence than CG for large systems)
- [ ] CHARMM22 force field (second FF, enables cross-validation)
- [ ] EEF1 implicit solvation (makes in-vacuo MD physically meaningful)
- [ ] MMFF94 force field (for small molecules/ligands)
- [ ] Nonbonded cutoff with neighbor lists (O(N) instead of O(N²))

### MD Extensions
- [ ] NPT ensemble (Berendsen or Parrinello-Rahman barostat)
- [ ] Constraint algorithms (SHAKE/LINCS for H-bond constraints, larger timestep)
- [ ] Trajectory output (DCD/XTC format)
- [ ] Replica exchange / parallel tempering

### Pipeline & Data
- [ ] `ferritin.prepare(structure)` convenience function (load → reconstruct → add H → minimize)
- [ ] Batch structure preparation pipeline (load_and_prepare)
- [ ] Per-residue feature extraction to Arrow (SASA, SS, dihedrals, H-bonds in one pass)
- [ ] Structural alphabet (3Di / VQ-VAE based)
- [ ] Interface residue detection (SASA in complex vs monomer)

### Performance
- [ ] Nonbonded neighbor list (cell list or Verlet list) for O(N) scaling
- [ ] SIMD vectorization for distance/energy computation
- [ ] GPU acceleration (CUDA/Metal) for nonbonded interactions — only when profiling proves need

### Python API
- [ ] `ferritin.prepare()` high-level pipeline function
- [ ] Trajectory analysis module (RMSD over time, Rg, etc.)
- [ ] Jupyter notebook examples
- [ ] PyPI release

---

## Design Principles

1. **Compute kernel, not platform** — zero friction on any infrastructure
2. **Batch-first** — every feature works on N structures with rayon parallelism
3. **Oracle-validated** — every feature tested against an independent C/C++ reference
4. **Pipeline-native** — Arrow/Parquet output, zero-GIL, composable
5. **No speculative abstractions** — build what's needed, test what's built

## Available Oracles

| Oracle | What | Status |
|--------|------|--------|
| FreeSASA (C) | SASA | 0.01% match with ProtOr |
| BALL / BiochemicalAlgorithms.jl | H placement, energy, reconstruction | Julia 1.11, working |
| Gemmi (C++) | Atom counts, H placement (needs full CCD) | CLI working |
| PDBFixer (OpenMM) | Full structure preparation | Python, working |
| GROMACS | DSSP, H-bonds, MD | Compiled, needs format conversion |
| C++ USalign | Alignment | Binary available |
