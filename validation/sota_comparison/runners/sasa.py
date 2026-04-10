"""SASA runners.

Per-op payload schema (`payload` field of RunnerResult):

    {
        "total_sasa": float,            # Total SASA in Å²
        "per_residue": [                # One entry per residue, ordered by
            {                           # (chain, resi, icode) — but JOIN ON
                "chain": "A",           # the (chain, resi, icode) tuple in
                "resi": 1,              # the aggregator, never positional.
                "icode": "",
                "name": "THR",          # 3-letter residue name
                "sasa": 12.34,          # Å²
            },
            ...
        ],
        "n_residues": int,
        "radii": "protor",              # which radii table was used
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
except ImportError as e:
    _FERRITIN_OK = False
    _FERRITIN_VERSION = ""
    _FERRITIN_IMPORT_ERROR = str(e)


if _FERRITIN_OK:

    @register("sasa", "ferritin")
    def ferritin(pdb_path: str) -> RunnerResult:
        """Compute SASA via ferritin.total_sasa + ferritin.residue_sasa.

        Uses radii="protor" to match FreeSASA's NACCESS-style table. Probe
        radius defaults to 1.4 Å (water), n_points defaults to 960 (matches
        ferritin's batch_total_sasa default).
        """
        s, _ = time_call(_ferritin.load, pdb_path)
        total, elapsed_total = time_call(
            _ferritin.total_sasa, s, radii="protor"
        )
        per_res_arr, elapsed_res = time_call(
            _ferritin.residue_sasa, s, radii="protor"
        )

        # Build the per-residue list keyed for downstream join.
        residues = list(s.residues)
        per_residue = []
        if len(residues) != len(per_res_arr):
            # Length mismatch — should not happen if ferritin's API is consistent.
            return RunnerResult(
                op="sasa",
                impl="ferritin",
                impl_version=_FERRITIN_VERSION,
                pdb_id="",
                pdb_path=pdb_path,
                elapsed_s=elapsed_total + elapsed_res,
                status="error",
                error=(
                    f"residue count mismatch: structure.residues has "
                    f"{len(residues)}, residue_sasa has {len(per_res_arr)}"
                ),
                payload={},
            )

        for r, val in zip(residues, per_res_arr):
            per_residue.append({
                "chain": r.chain_id,
                "resi": int(r.serial_number),
                "icode": r.insertion_code or "",
                "name": r.name or "",
                "sasa": float(val),
            })

        return RunnerResult(
            op="sasa",
            impl="ferritin",
            impl_version=_FERRITIN_VERSION,
            pdb_id="",  # filled in by the driver
            pdb_path=pdb_path,
            elapsed_s=elapsed_total + elapsed_res,
            status="ok",
            error=None,
            payload={
                "total_sasa": float(total),
                "per_residue": per_residue,
                "n_residues": len(per_residue),
                "radii": "protor",
            },
        )


# ---------------------------------------------------------------------------
# FreeSASA
# ---------------------------------------------------------------------------

try:
    import freesasa as _freesasa  # noqa: F401  (only used inside the runner)
    _FREESASA_OK = True
    try:
        _FREESASA_VERSION = "freesasa " + _metadata.version("freesasa")
    except Exception:
        _FREESASA_VERSION = "freesasa (unknown version)"
except ImportError:
    _FREESASA_OK = False
    _FREESASA_VERSION = ""


if _FREESASA_OK:

    @register("sasa", "freesasa")
    def freesasa(pdb_path: str) -> RunnerResult:
        """Compute SASA via the freesasa Python package.

        Matches ferritin's defaults as closely as possible:

        - Include HETATM atoms (ferritin loads them; freesasa default skips).
          This is the single biggest contributor to the agreement: 1pgb alone
          has 24 waters that ferritin counts and stock freesasa doesn't.
        - Use the ProtOr classifier explicitly (FreeSASA default is OONS-like
          Lee-Richards; ferritin's `radii="protor"` should match ProtOr).
        - Use the Shrake-Rupley algorithm with 960 points (matches ferritin's
          default).
        - Probe radius 1.4 Å (matches ferritin default).

        Per-residue results are keyed by (chain, resi, icode) tuples and
        normalized to the same shape as the ferritin runner so the aggregator
        can join them directly.
        """
        options = {
            "hetatm": True,          # include HETATM to match ferritin
            "hydrogen": False,       # ferritin-loaded crystal files have no H
            "join-models": False,    # first model only, same as ferritin
            "skip-unknown": False,
            "halt-at-unknown": False,
        }
        classifier = _freesasa.Classifier.getStandardClassifier("protor")
        structure, elapsed_load = time_call(
            _freesasa.Structure, pdb_path, classifier, options
        )
        params = _freesasa.Parameters(
            {
                "algorithm": _freesasa.ShrakeRupley,
                "n-points": 960,
                "probe-radius": 1.4,
            }
        )
        result, elapsed_calc = time_call(_freesasa.calc, structure, params)
        total = result.totalArea()

        # Per-residue: residueAreas() returns dict[chain_id -> dict[res_num_str -> ResidueArea]]
        # The res_num_str may include insertion code as e.g. "65A".
        per_residue = []
        residue_areas = result.residueAreas()
        for chain_id, by_resnum in residue_areas.items():
            for resnum_str, ra in by_resnum.items():
                # Split off insertion code if present.
                icode = ""
                resi_str = resnum_str
                if resi_str and not (resi_str[-1].isdigit() or resi_str[-1] == "-"):
                    icode = resi_str[-1]
                    resi_str = resi_str[:-1]
                try:
                    resi = int(resi_str)
                except ValueError:
                    # Can't parse — skip but record.
                    continue
                per_residue.append({
                    "chain": chain_id,
                    "resi": resi,
                    "icode": icode,
                    "name": ra.residueType,
                    "sasa": float(ra.total),
                })

        return RunnerResult(
            op="sasa",
            impl="freesasa",
            impl_version=_FREESASA_VERSION,
            pdb_id="",
            pdb_path=pdb_path,
            elapsed_s=elapsed_load + elapsed_calc,
            status="ok",
            error=None,
            payload={
                "total_sasa": float(total),
                "per_residue": per_residue,
                "n_residues": len(per_residue),
                "radii": "freesasa-default",
            },
        )
