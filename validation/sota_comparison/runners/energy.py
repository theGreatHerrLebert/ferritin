"""Energy runners.

Per-op payload schema (`payload` field of RunnerResult):

    {
        "total": float,                # Total energy in kJ/mol
        "units": "kJ/mol",
        "ff": "amber96",               # which force field was used
        "components": {                # Per-component energy in kJ/mol;
            "bond_stretch": float,     # None for components the FF doesn't
            "angle_bend": float,       # define (e.g. solvation in vacuum FFs).
            "torsion": float,
            "improper_torsion": float,
            "vdw": float,
            "electrostatic": float,
            "solvation": Optional[float],
        },
        "n_unassigned_atoms": int,     # ferritin reports this; -1 for tools
                                       # that don't.
    }
"""

from __future__ import annotations

import importlib.metadata as _metadata
from typing import Optional

from ._base import (
    RunnerResult,
    register,
    time_call,
)


# ---------------------------------------------------------------------------
# ferritin baseline
# ---------------------------------------------------------------------------

try:
    import ferritin as _ferritin
    _FERRITIN_OK = True
    try:
        _FERRITIN_VERSION = "ferritin " + _metadata.version("ferritin")
    except Exception:
        _FERRITIN_VERSION = "ferritin (unknown version)"
except ImportError:
    _FERRITIN_OK = False
    _FERRITIN_VERSION = ""


_COMPONENT_KEYS = (
    "bond_stretch",
    "angle_bend",
    "torsion",
    "improper_torsion",
    "vdw",
    "electrostatic",
    "solvation",
)


def _normalize_components(d: dict) -> dict:
    """Pick out the canonical component keys from a ferritin energy dict.

    Returns a dict with every key from _COMPONENT_KEYS present, using None for
    components the source dict doesn't have. Defensive against ferritin adding
    new components in the future without breaking the schema.
    """
    return {k: float(d[k]) if d.get(k) is not None else None for k in _COMPONENT_KEYS}


if _FERRITIN_OK:

    @register("energy", "ferritin")
    def ferritin(pdb_path: str) -> RunnerResult:
        """Compute AMBER96 energy via ferritin.compute_energy.

        Returns kJ/mol (ferritin's default unit). Used as the baseline against
        which OpenMM and BALL are compared.
        """
        s, _ = time_call(_ferritin.load, pdb_path)
        result, elapsed = time_call(
            _ferritin.compute_energy, s, ff="amber96", units="kJ/mol"
        )

        return RunnerResult(
            op="energy",
            impl="ferritin",
            impl_version=_FERRITIN_VERSION,
            pdb_id="",
            pdb_path=pdb_path,
            elapsed_s=elapsed,
            status="ok",
            error=None,
            payload={
                "total": float(result["total"]),
                "units": "kJ/mol",
                "ff": "amber96",
                "components": _normalize_components(result),
                "n_unassigned_atoms": int(result.get("n_unassigned_atoms", 0)),
            },
        )


# ---------------------------------------------------------------------------
# OpenMM (added in Step 4)
# ---------------------------------------------------------------------------
# Placeholder: the openmm runner is wired in Step 4 of the implementation
# plan, after the FreeSASA / SASA pipeline has shipped end-to-end.
