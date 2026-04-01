"""Pythonic wrappers for structure I/O."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .core import RustWrapperObject

try:
    import ferritin_connector
    ims = ferritin_connector.py_structure
except ImportError:
    ims = None


class StructureData(RustWrapperObject):
    """Loaded structure data from a PDB or mmCIF file.

    Attributes:
        coords: CA/C3' coordinates as Nx3 numpy array.
        sequence: One-letter residue sequence.
        sec_structure: Secondary structure string (H/E/T/C).
        chain_id: Chain identifier.
    """

    def __init__(self, ptr):
        self._ptr = ptr

    @classmethod
    def from_py_ptr(cls, ptr) -> StructureData:
        return cls(ptr)

    def get_py_ptr(self):
        return self._ptr

    @property
    def coords(self) -> NDArray[np.float64]:
        return self._ptr.coords

    @property
    def sequence(self) -> str:
        return self._ptr.sequence

    @property
    def sec_structure(self) -> str:
        return self._ptr.sec_structure

    @property
    def chain_id(self) -> str:
        return self._ptr.chain_id

    def __len__(self) -> int:
        return len(self._ptr)

    def __repr__(self) -> str:
        return f"StructureData(chain='{self.chain_id}', n_residues={len(self)})"
