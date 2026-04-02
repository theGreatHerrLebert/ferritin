"""DSSP secondary structure assignment (Kabsch-Sander algorithm).

Native Rust implementation — no external DSSP binary needed.

Assigns per-residue secondary structure from backbone coordinates using
the Kabsch-Sander hydrogen bond energy criterion.

Functions:
    dssp                — SS string (H/G/I/E/B/T/S/C per residue)
    dssp_array          — SS as numpy u8 array
    batch_dssp          — SS for many structures in parallel
    load_and_dssp       — load files + compute SS in one parallel call
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

try:
    import ferritin_connector
    _dssp = ferritin_connector.py_dssp
except ImportError:
    _dssp = None


def _get_ptr(structure):
    if hasattr(structure, 'get_py_ptr'):
        return structure.get_py_ptr()
    return structure


def dssp(structure) -> str:
    """Assign DSSP secondary structure.

    Uses the Kabsch-Sander hydrogen bond energy criterion to assign
    one-letter SS codes per amino acid residue:
        H = alpha helix (4-turn helix)
        G = 3-10 helix (3-turn helix)
        I = pi helix (5-turn helix)
        E = extended strand (part of beta ladder)
        B = isolated beta bridge
        T = hydrogen-bonded turn
        S = bend (CA angle > 70 degrees)
        C = coil (none of the above)

    No external DSSP binary needed — pure Rust implementation.

    Args:
        structure: A ferritin Structure.

    Returns:
        String of SS codes, one per amino acid residue.

    Examples:
        >>> ss = ferritin.dssp(structure)
        >>> print(ss)
        'CCCSSHHHHHHHHHHHCCCTHHHHHHTCCEEEEECCCCCCCTEEEEC'
        >>> n_helix = ss.count('H')

    Agent Notes:
        OUTPUT: String of length = number of amino acid residues (not total
        residues — water/ligands are excluded). Length may differ from
        structure.residue_count if non-AA residues are present.

        INTERPRET: H=alpha helix, E=beta strand are the main secondary
        structure types. G=3-10 helix and I=pi helix are rare but valid.
        B=isolated beta bridge (not part of a sheet). T=turn, S=bend, C=coil.

        REQUIRES: Backbone atoms N, CA, C, O must be present per residue.
        Residues missing any of these are skipped. Structures with only CA
        atoms (e.g., coarse-grained models) will return empty strings —
        use assign_secondary_structure() instead for CA-only approximation.

        COMPARE: Results may differ slightly from the DSSP binary (mkdssp)
        due to implementation differences in bridge detection and edge cases.
        Core helix/strand regions match well; disagreements are typically
        at boundaries (1-2 residues).

        PREFER: For many structures, use batch_dssp(structures, n_threads=-1).
    """
    return _dssp.dssp(_get_ptr(structure))


def dssp_array(structure) -> NDArray[np.uint8]:
    """Assign DSSP as a numpy array of ASCII character codes.

    Same as dssp() but returns numpy array for vectorized operations.

    Examples:
        >>> ss = ferritin.dssp_array(structure)
        >>> is_helix = ss == ord('H')
        >>> helix_fraction = is_helix.mean()
    """
    return np.asarray(_dssp.dssp_array(_get_ptr(structure)))


def batch_dssp(
    structures: Sequence,
    *,
    n_threads: Optional[int] = None,
) -> List[str]:
    """Compute DSSP for many structures in parallel (Rust + rayon).

    Args:
        structures: List of ferritin Structure objects.
        n_threads: Thread count. None/-1 = all cores.

    Returns:
        List of SS strings.
    """
    ptrs = [_get_ptr(s) for s in structures]
    return _dssp.batch_dssp(ptrs, n_threads)


def load_and_dssp(
    paths: Sequence,
    *,
    n_threads: Optional[int] = None,
) -> List[Tuple[int, str]]:
    """Load files and compute DSSP in one parallel call (zero GIL).

    Args:
        paths: List of file paths.
        n_threads: Thread count. None/-1 = all cores.

    Returns:
        List of (index, ss_string) tuples.
    """
    str_paths = [str(p) for p in paths]
    return _dssp.load_and_dssp(str_paths, n_threads)
