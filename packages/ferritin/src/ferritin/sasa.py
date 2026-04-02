"""Solvent Accessible Surface Area (SASA) — Shrake-Rupley algorithm.

Rust implementation with rayon parallelism. Computes per-atom and per-residue
SASA using the Shrake-Rupley numerical dot method.

Functions:
    atom_sasa          — per-atom SASA (Angstroms²)
    residue_sasa       — per-residue SASA (sum of atom contributions)
    relative_sasa      — RSA (residue SASA / max SASA for residue type)
    total_sasa         — total SASA of a structure
    batch_total_sasa   — total SASA for many structures in parallel
    load_and_sasa      — load files + compute SASA in one parallel call
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

try:
    import ferritin_connector

    _sasa = ferritin_connector.py_sasa
except ImportError:
    _sasa = None


def _get_ptr(structure):
    if hasattr(structure, 'get_py_ptr'):
        return structure.get_py_ptr()
    return structure


def atom_sasa(
    structure,
    probe: float = 1.4,
    n_points: int = 960,
) -> NDArray[np.float64]:
    """Compute per-atom SASA using Shrake-Rupley algorithm.

    Args:
        structure: A ferritin Structure.
        probe: Probe radius in Angstroms (default 1.4 for water).
        n_points: Test points per sphere (default 960; higher = more precise).

    Returns:
        1D numpy array of per-atom SASA in Angstroms².

    Examples:
        >>> sasa = ferritin.atom_sasa(structure)
        >>> print(f"Total SASA: {sasa.sum():.0f} A²")
        >>> exposed = sasa > 0
        >>> print(f"{exposed.sum()} exposed atoms out of {len(sasa)}")

    Agent Notes:
        DEFAULTS: probe=1.4 is the standard water probe radius. Don't change
        unless you have a specific reason (e.g., probe=0 for van der Waals surface).

        PRECISION: n_points=960 gives ~0.2% accuracy. Use 100 for fast screening,
        960 for publication quality. Doubling points halves the error but doubles time.

        OUTPUT: Array length = structure.atom_count (all atoms including HETATM/water).
        To get per-residue or per-chain SASA, use residue_sasa() instead.

        PREFER: For many structures, use batch_total_sasa() with n_threads=-1.

        COST: O(N * P * k) where N=atoms, P=points, k=avg neighbors.
        Crambin (327 atoms): ~12ms. Large complex (58k atoms): ~230ms.
    """
    return np.asarray(_sasa.atom_sasa(_get_ptr(structure), probe, n_points))


def residue_sasa(
    structure,
    probe: float = 1.4,
    n_points: int = 960,
) -> NDArray[np.float64]:
    """Compute per-residue SASA (sum of atom contributions per residue).

    Args:
        structure: A ferritin Structure.
        probe: Probe radius in Angstroms (default 1.4).
        n_points: Test points per sphere (default 960).

    Returns:
        1D numpy array of per-residue SASA in Angstroms².
    """
    return np.asarray(_sasa.residue_sasa(_get_ptr(structure), probe, n_points))


def relative_sasa(
    structure,
    probe: float = 1.4,
    n_points: int = 960,
) -> NDArray[np.float64]:
    """Compute relative solvent accessibility (RSA) per residue.

    RSA = residue_SASA / max_SASA_for_residue_type (Tien et al. 2013).
    Values > 0.25 are typically considered "exposed".

    Args:
        structure: A ferritin Structure.
        probe: Probe radius in Angstroms (default 1.4).
        n_points: Test points per sphere (default 960).

    Returns:
        1D numpy array of RSA values (0.0–1.0+). NaN for non-standard residues.

    Agent Notes:
        INTERPRET: RSA < 0.25 = buried (core), RSA >= 0.25 = exposed (surface).
        This threshold is the standard in literature (Tien et al. 2013).

        WATCH: RSA can exceed 1.0 for residues in extended conformations or
        at chain termini. This is normal, not an error.

        WATCH: NaN values indicate non-standard residues (ligands, modified
        amino acids, water) for which no reference max SASA exists.
        Filter with ~np.isnan(rsa) before analysis.

        USE FOR: Burial classification, identifying surface residues for
        mutation studies, interface residue detection (compare RSA in
        complex vs monomer).
    """
    return np.asarray(_sasa.relative_sasa(_get_ptr(structure), probe, n_points))


def total_sasa(
    structure,
    probe: float = 1.4,
    n_points: int = 960,
) -> float:
    """Total SASA of a structure in Angstroms².

    Args:
        structure: A ferritin Structure.
        probe: Probe radius in Angstroms (default 1.4).
        n_points: Test points per sphere (default 960).

    Returns:
        Total SASA in Angstroms².

    Examples:
        >>> print(f"Total SASA: {ferritin.total_sasa(structure):.0f} A²")
    """
    return _sasa.total_sasa(_get_ptr(structure), probe, n_points)


def batch_total_sasa(
    structures: Sequence,
    probe: float = 1.4,
    n_points: int = 960,
    *,
    n_threads: Optional[int] = None,
) -> NDArray[np.float64]:
    """Compute total SASA for many structures in parallel (Rust + rayon).

    Args:
        structures: List of ferritin Structure objects.
        probe: Probe radius in Angstroms (default 1.4).
        n_points: Test points per sphere (default 960).
        n_threads: Thread count. None/-1 = all cores.

    Returns:
        1D numpy array of total SASA values.
    """
    ptrs = [_get_ptr(s) for s in structures]
    return np.asarray(_sasa.batch_total_sasa(ptrs, probe, n_points, n_threads))


def load_and_sasa(
    paths: Sequence,
    probe: float = 1.4,
    n_points: int = 960,
    *,
    n_threads: Optional[int] = None,
) -> List[Tuple[int, float]]:
    """Load files and compute total SASA in one parallel call (zero GIL).

    Args:
        paths: List of file paths.
        probe: Probe radius (default 1.4).
        n_points: Test points per sphere (default 960).
        n_threads: Thread count. None/-1 = all cores.

    Returns:
        List of (index, total_sasa) tuples for files that loaded.

    Examples:
        >>> results = ferritin.load_and_sasa(pdb_files, n_threads=-1)
        >>> for idx, sasa in results:
        ...     print(f"{pdb_files[idx]}: {sasa:.0f} A²")
    """
    str_paths = [str(p) for p in paths]
    return _sasa.load_and_sasa(str_paths, probe, n_points, n_threads)
