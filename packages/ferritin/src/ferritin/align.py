"""Pythonic wrappers for structural alignment."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .core import RustWrapperObject

try:
    import ferritin_connector
    ims = ferritin_connector.py_align
except ImportError:
    ims = None


class AlignResult(RustWrapperObject):
    """Result of a structural alignment.

    Attributes:
        tm_score_chain1: TM-score normalized by length of chain 1.
        tm_score_chain2: TM-score normalized by length of chain 2.
        rmsd: RMSD of aligned residues.
        n_aligned: Number of aligned residue pairs.
        seq_identity: Fraction of identical aligned residues.
        rotation_matrix: 3x3 rotation matrix (numpy array).
        translation: Translation vector [x, y, z] (numpy array).
        aligned_seq_x: Aligned sequence of structure 1 (with gaps).
        aligned_seq_y: Aligned sequence of structure 2 (with gaps).
    """

    def __init__(self, ptr):
        self._ptr = ptr

    @classmethod
    def from_py_ptr(cls, ptr) -> AlignResult:
        return cls(ptr)

    def get_py_ptr(self):
        return self._ptr

    @property
    def tm_score_chain1(self) -> float:
        return self._ptr.tm_score_chain1

    @property
    def tm_score_chain2(self) -> float:
        return self._ptr.tm_score_chain2

    @property
    def rmsd(self) -> float:
        return self._ptr.rmsd

    @property
    def n_aligned(self) -> int:
        return self._ptr.n_aligned

    @property
    def seq_identity(self) -> float:
        return self._ptr.seq_identity

    @property
    def rotation_matrix(self) -> NDArray[np.float64]:
        return self._ptr.rotation_matrix

    @property
    def translation(self) -> NDArray[np.float64]:
        return np.array(self._ptr.translation)

    @property
    def aligned_seq_x(self) -> str:
        return self._ptr.aligned_seq_x

    @property
    def aligned_seq_y(self) -> str:
        return self._ptr.aligned_seq_y

    def __repr__(self) -> str:
        return (
            f"AlignResult(TM1={self.tm_score_chain1:.4f}, "
            f"TM2={self.tm_score_chain2:.4f}, "
            f"RMSD={self.rmsd:.2f}, "
            f"Lali={self.n_aligned})"
        )
