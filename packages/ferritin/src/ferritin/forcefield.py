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
    """
    return _ff.minimize_hydrogens(_get_ptr(structure), max_steps, gradient_tolerance)


def minimize_structure(
    structure,
    max_steps: int = 1000,
    gradient_tolerance: float = 0.1,
) -> dict:
    """Full energy minimization using AMBER96 force field.

    All atoms are free to move.

    Args:
        structure: A ferritin Structure.
        max_steps: Maximum optimization steps (default 1000).
        gradient_tolerance: Convergence criterion in kcal/mol/A (default 0.1).

    Returns:
        dict with: coords (Nx3), initial_energy, final_energy,
        steps, converged, energy_components.
    """
    return _ff.minimize_structure(_get_ptr(structure), max_steps, gradient_tolerance)
