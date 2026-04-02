"""AMBER force field energy computation and minimization.

Implements AMBER96 energy terms (bond stretch, angle bend, torsion,
Lennard-Jones, Coulomb) with steepest descent minimization.

Functions:
    compute_energy      — AMBER energy breakdown
    minimize_hydrogens  — optimize H positions (freeze heavy atoms)
    minimize_structure  — full energy minimization
"""

from __future__ import annotations

try:
    import ferritin_connector
    _ff = ferritin_connector.py_forcefield
except ImportError:
    _ff = None


def _get_ptr(structure):
    if hasattr(structure, 'get_py_ptr'):
        return structure.get_py_ptr()
    return structure


def compute_energy(structure) -> dict:
    """Compute AMBER96 force field energy of a structure.

    Returns dict with energy components (all in kcal/mol):
        bond_stretch, angle_bend, torsion, vdw, electrostatic, total

    Examples:
        >>> e = ferritin.compute_energy(structure)
        >>> print(f"Total: {e['total']:.1f} kcal/mol")

    Agent Notes:
        INTERPRET: Absolute energy values are NOT meaningful across structures
        of different sizes. Only compare energies of the SAME structure before/after
        modification (e.g., before/after minimization).

        WATCH: vdw > 1e6 usually means steric clashes (atoms overlapping).
        This is common in raw X-ray structures before minimization — it is NOT a bug.

        WATCH: electrostatic energy is large without solvent. In vacuo energies
        are always higher than solvated. Don't interpret large electrostatic
        values as errors.

        EXPECT: For a well-refined crystal structure like crambin (46 residues):
        bond_stretch ~80, angle_bend ~70, torsion ~140, vdw ~2400,
        electrostatic ~20000 kcal/mol.

        COST: O(N^2) for nonbonded terms. Structures > 5000 atoms take seconds.
    """
    return _ff.compute_energy(_get_ptr(structure))


def minimize_hydrogens(
    structure,
    max_steps: int = 500,
    gradient_tolerance: float = 0.1,
) -> dict:
    """Minimize hydrogen positions using AMBER96 force field.

    Freezes all heavy atoms and optimizes only H positions.
    Uses steepest descent with adaptive step size.

    Args:
        structure: A ferritin Structure.
        max_steps: Maximum optimization steps (default 500).
        gradient_tolerance: Convergence criterion in kcal/mol/A (default 0.1).

    Returns:
        dict with: coords (Nx3), initial_energy, final_energy,
        steps, converged, energy_components.

    Examples:
        >>> result = ferritin.minimize_hydrogens(structure)
        >>> print(f"Energy: {result['initial_energy']:.0f} -> {result['final_energy']:.0f}")

    Agent Notes:
        PREREQUISITE: Structure must already contain hydrogen atoms. Most PDB
        files from X-ray crystallography do NOT have hydrogens. If atom_count
        seems low or there are no H atoms, hydrogens need to be added first.

        VERIFY: Always check result['converged']. If False, the minimization
        hit max_steps without converging — increase max_steps or check for
        severe clashes.

        EXPECT: For refined crystal structures, this often converges in 1 step
        (hydrogens are already well-placed). Large energy changes indicate the
        input had bad hydrogen positions.

        PREFER: Use batch_minimize_hydrogens() when processing multiple structures.
        It uses rayon parallelism across structures.

        COST: O(N^2) per step for nonbonded interactions. Fast for small proteins
        (<1000 atoms), slower for large complexes.
    """
    return _ff.minimize_hydrogens(_get_ptr(structure), max_steps, gradient_tolerance)


def minimize_structure(
    structure,
    max_steps: int = 1000,
    gradient_tolerance: float = 0.1,
) -> dict:
    """Full energy minimization using AMBER96 force field.

    All atoms are free to move. Uses steepest descent.

    Args:
        structure: A ferritin Structure.
        max_steps: Maximum optimization steps (default 1000).
        gradient_tolerance: Convergence criterion in kcal/mol/A (default 0.1).

    Returns:
        dict with: coords (Nx3), initial_energy, final_energy,
        steps, converged, energy_components.

    Agent Notes:
        CAUTION: Full minimization moves ALL atoms. This can distort the
        structure if run for too many steps. For most use cases,
        minimize_hydrogens() is safer and sufficient.

        CAUTION: Steepest descent is robust but slow to converge. It will
        resolve major clashes but won't find the true energy minimum.
        This is appropriate for structure preparation, NOT for finding
        native conformations.

        WATCH: If the structure deforms (RMSD > 2A from input), reduce
        max_steps. Typical use: 100-500 steps to relieve clashes.

        AVOID: Don't use this to "refine" a docked pose or homology model
        beyond removing clashes. The in-vacuo force field will collapse
        exposed sidechains without solvent.
    """
    return _ff.minimize_structure(_get_ptr(structure), max_steps, gradient_tolerance)


def batch_minimize_hydrogens(
    structures,
    max_steps: int = 500,
    gradient_tolerance: float = 0.1,
    *,
    n_threads=None,
):
    """Minimize hydrogen positions for many structures in parallel (Rust + rayon).

    Args:
        structures: List of ferritin Structure objects.
        max_steps: Maximum optimization steps per structure.
        gradient_tolerance: Convergence criterion in kcal/mol/A.
        n_threads: Thread count. None/-1 = all cores.

    Returns:
        List of result dicts (same format as minimize_hydrogens).

    Examples:
        >>> results = ferritin.batch_minimize_hydrogens(structures, n_threads=-1)
        >>> for r in results:
        ...     print(f"E: {r['initial_energy']:.0f} -> {r['final_energy']:.0f}")

    Agent Notes:
        PREFER: This over a Python loop — uses rayon for true multi-core
        parallelism with GIL released.

        SCALE: Speedup scales with number of structures and cores.
        Best for datasets of 10+ structures on multi-core machines.
    """
    ptrs = [_get_ptr(s) for s in structures]
    return _ff.batch_minimize_hydrogens(ptrs, max_steps, gradient_tolerance, n_threads)


def load_and_minimize_hydrogens(
    paths,
    max_steps: int = 500,
    gradient_tolerance: float = 0.1,
    *,
    n_threads=None,
):
    """Load files and minimize hydrogens in one parallel call (zero GIL).

    Entire pipeline (I/O, parsing, topology, minimization) runs in Rust.

    Args:
        paths: List of file paths.
        max_steps: Maximum optimization steps per structure.
        gradient_tolerance: Convergence criterion in kcal/mol/A.
        n_threads: Thread count. None/-1 = all cores.

    Returns:
        List of (index, result_dict) tuples. Files that fail to load are skipped.

    Agent Notes:
        PREFER: This is the fastest way to minimize many structures. Zero Python
        overhead — everything happens in Rust with rayon parallelism.

        TOLERANT: Files that fail to load are silently skipped. Check the
        returned indices to see which files succeeded.
    """
    str_paths = [str(p) for p in paths]
    return _ff.load_and_minimize_hydrogens(str_paths, max_steps, gradient_tolerance, n_threads)
