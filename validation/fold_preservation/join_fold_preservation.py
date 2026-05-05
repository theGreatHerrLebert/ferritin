"""Join proteon and OpenMM fold-preservation runs into per-pair JSONLs.

Each input run is a 1k-PDB sweep where one minimizer (proteon CHARMM19+EEF1
or proteon AMBER96 or OpenMM CHARMM36+OBC2 or OpenMM AMBER96) ran on a
seeded random sample. The runs share the same input corpus + seed so we
can join records by ``pdb`` to compare TM-score / RMSD between the two
implementations on each structure.

Outputs two joined JSONLs:

    validation/fold_preservation/charmm_pair_1k.jsonl
        proteon CHARMM19+EEF1  vs  OpenMM CHARMM36+OBC2

    validation/fold_preservation/amber_pair_1k.jsonl
        proteon AMBER96  vs  OpenMM AMBER96

Each line:

    {
      "pdb": "1abc.pdb",
      "n_ca": 142,
      "n_ca_pre": 142,
      "n_ca_post_proteon": 142,
      "n_ca_post_openmm": 142,
      "proteon": {"tm_score": 0.9943, "rmsd": 0.31,
                  "initial_energy": ..., "final_energy": ...,
                  "minimizer_steps": 100, "converged": false},
      "openmm":  {"tm_score": 0.9991, "rmsd": 0.18,
                  "initial_energy_kj": ..., "final_energy_kj": ...,
                  "wall_s": 28.3},
      "tm_diff": 0.0048,            // openmm - proteon (positive = OpenMM holds the fold tighter)
      "rmsd_diff_A": -0.13           // openmm - proteon
    }

Skipped / errored records on EITHER side mark the joined record as
{"pdb": ..., "skipped": "<which-side>: <reason>"}.

Reproducing this artifact:

    1. Run proteon CHARMM19+EEF1 on the 1k seeded sample ⇒ tm_fold_preservation.jsonl
    2. Run OpenMM CHARMM36+OBC2 on the same sample      ⇒ tm_fold_preservation_openmm.jsonl
    3. python validation/fold_preservation/join_fold_preservation.py

The first two steps are heavy (hours of GPU + CPU compute on monster3).
The join itself is O(N) and seconds.
"""
from __future__ import annotations

import json
import os
import pathlib
import sys


HERE = pathlib.Path(__file__).resolve().parent

# v0.2.0 data-mount contract: when PROTEON_OUTPUT_DIR is set the joiner reads
# the per-side JSONLs from that directory (matching where the runners wrote
# them) and emits the joined artifacts to the same place. When unset, fall
# back to the historical HERE-based layout used by the monster3-→-repo scp
# workflow.
_OUT = os.environ.get("PROTEON_OUTPUT_DIR")
SRC_DIR = pathlib.Path(_OUT) if _OUT else HERE
DEST_DIR = pathlib.Path(_OUT) if _OUT else HERE


def _load_by_pdb(path: pathlib.Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            d = json.loads(line)
            pdb = d.get("pdb")
            if pdb:
                out[pdb] = d
    return out


def _classify(rec: dict) -> str:
    if not rec:
        return "missing"
    if rec.get("error"):
        return "error"
    if rec.get("skipped"):
        return "skipped"
    return "ok"


def _join_pair(
    proteon_jsonl: pathlib.Path,
    openmm_jsonl: pathlib.Path,
    out_jsonl: pathlib.Path,
    proteon_label: str,
    openmm_label: str,
) -> None:
    """Join the two side-runs by pdb, emit one canonical record per PDB."""
    p = _load_by_pdb(proteon_jsonl)
    o = _load_by_pdb(openmm_jsonl)
    pdbs = sorted(set(p.keys()) | set(o.keys()))

    n_ok = n_skip = n_err = 0
    out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with out_jsonl.open("w", encoding="utf-8") as f:
        for pdb in pdbs:
            pr = p.get(pdb) or {}
            om = o.get(pdb) or {}
            cp = _classify(pr)
            co = _classify(om)

            if cp != "ok" or co != "ok":
                rec = {
                    "pdb": pdb,
                    "skipped": f"proteon={cp}; openmm={co}",
                }
                if cp == "error":
                    rec["proteon_error"] = pr.get("error")
                if co == "error":
                    rec["openmm_error"] = om.get("error")
                n_skip += 1 if (cp != "error" and co != "error") else 0
                n_err += 1 if (cp == "error" or co == "error") else 0
                f.write(json.dumps(rec) + "\n")
                continue

            n_ok += 1
            tm_p = float(pr["tm_score"])
            tm_o = float(om["tm_score"])
            rmsd_p = float(pr["rmsd"])
            rmsd_o = float(om["rmsd"])
            rec = {
                "pdb": pdb,
                "n_ca": int(pr.get("n_ca") or om.get("n_ca") or 0),
                "n_ca_pre": int(pr.get("n_ca_pre") or om.get("n_ca_pre") or 0),
                "n_ca_post_proteon": int(pr.get("n_ca_post") or 0),
                "n_ca_post_openmm": int(om.get("n_ca_post") or 0),
                "proteon": {
                    "tm_score": tm_p,
                    "rmsd": rmsd_p,
                    "initial_energy": pr.get("initial_energy"),
                    "final_energy": pr.get("final_energy"),
                    "minimizer_steps": pr.get("minimizer_steps"),
                    "converged": pr.get("converged"),
                    "label": proteon_label,
                },
                "openmm": {
                    "tm_score": tm_o,
                    "rmsd": rmsd_o,
                    "initial_energy_kj": om.get("initial_energy_kj"),
                    "final_energy_kj": om.get("final_energy_kj"),
                    "wall_s": om.get("wall_s"),
                    "label": openmm_label,
                },
                "tm_diff": tm_o - tm_p,
                "rmsd_diff_A": rmsd_o - rmsd_p,
            }
            f.write(json.dumps(rec) + "\n")

    print(
        f"{out_jsonl.name}: ok={n_ok} skip={n_skip} err={n_err} "
        f"of {len(pdbs)} PDBs"
    )


def main() -> int:
    pairs = [
        (
            SRC_DIR / "tm_fold_preservation.jsonl",
            SRC_DIR / "tm_fold_preservation_openmm.jsonl",
            DEST_DIR / "charmm_pair_1k.jsonl",
            "proteon CHARMM19+EEF1",
            "OpenMM CHARMM36+OBC2",
        ),
        (
            SRC_DIR / "tm_fold_preservation_amber.jsonl",
            SRC_DIR / "tm_fold_preservation_openmm_amber.jsonl",
            DEST_DIR / "amber_pair_1k.jsonl",
            "proteon AMBER96",
            "OpenMM AMBER96",
        ),
    ]
    missing = []
    for src1, src2, _, _, _ in pairs:
        for src in (src1, src2):
            if not src.is_file():
                missing.append(str(src))
    if missing:
        print(f"missing input JSONLs: {missing}", file=sys.stderr)
        return 2

    for src1, src2, out, lbl1, lbl2 in pairs:
        _join_pair(src1, src2, out, lbl1, lbl2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
