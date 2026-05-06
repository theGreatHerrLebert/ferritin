"""TM-align fold preservation benchmark — proteon AMBER96 variant.

Mirrors tm_fold_preservation.py but runs proteon with ff='amber96'
instead of charmm19_eef1. Produces a second dataset that lets us
compare proteon's two supported force fields side-by-side on the same
1000-PDB sample, and pairs with tm_fold_preservation_openmm_amber.py
for a three-tool AMBER96 comparison.
"""
from __future__ import annotations

import json
import os
import random
import time
import warnings
from pathlib import Path

import numpy as np

import proteon

# v0.2.0 data-mount contract: PROTEON_CORPUS_DIR / PROTEON_OUTPUT_DIR override
# the monster3 defaults. Set by the EVIDENT image entrypoint when the user
# bind-mounts /data/pdbs and /data/out.
PDB_DIR = Path(os.environ.get("PROTEON_CORPUS_DIR") or "/globalscratch/dateschn/proteon-benchmark/pdbs_50k")
OUTPUT_DIR = Path(os.environ.get("PROTEON_OUTPUT_DIR") or "/globalscratch/dateschn/proteon-benchmark")
OUT = OUTPUT_DIR / "tm_fold_preservation_amber.jsonl"
N = int(os.environ.get("N_PDBS", "1000"))
SEED = 42
CHUNK = 25
MINIMIZE_STEPS = 100

# Silence the AMBER96-experimental warning for batch runs.
warnings.filterwarnings("ignore", category=UserWarning, module="proteon.forcefield")


def tm_pair(ca_ref: np.ndarray, ca_mov: np.ndarray) -> dict:
    n = len(ca_ref)
    invmap = np.arange(n, dtype=np.int32)
    tm, n_aln, rmsd_val, _R, _t = proteon.tm_score(ca_mov, ca_ref, invmap)
    return {
        "tm_score": float(tm),
        "rmsd": float(rmsd_val),
        "n_ca": int(n),
        "n_aligned": int(n_aln),
    }


def process_chunk(paths_chunk, out_fh):
    results = []
    ref_loaded = proteon.batch_load_tolerant([str(p) for p in paths_chunk])
    ref_by_idx = dict(ref_loaded)
    ca_pre_by_idx = {i: proteon.extract_ca_coords(s) for i, s in ref_loaded}

    work_loaded = proteon.batch_load_tolerant([str(p) for p in paths_chunk])
    work_by_idx = dict(work_loaded)

    both_idx = sorted(set(ref_by_idx) & set(work_by_idx))
    missing_idx = sorted((set(range(len(paths_chunk))) - set(ref_by_idx))
                         | (set(range(len(paths_chunk))) - set(work_by_idx)))
    for mi in missing_idx:
        rec = {"pdb": paths_chunk[mi].name, "error": "load_failed"}
        out_fh.write(json.dumps(rec) + "\n")
        results.append(rec)

    if not both_idx:
        out_fh.flush()
        return results

    work_structs = [work_by_idx[i] for i in both_idx]
    try:
        # constrain_heavy=False: override AMBER96's H-only default so that
        # heavy atoms are free to move. Needed for apples-to-apples fold
        # preservation comparison with OpenMM/GROMACS (which minimize all
        # atoms). The FF-aware default is True for AMBER96 because heavy-
        # atom vacuum minimization is normally ill-posed; here we're
        # explicitly asking for it in the benchmark.
        reports = proteon.batch_prepare(
            work_structs, ff="amber96", minimize_steps=MINIMIZE_STEPS,
            constrain_heavy=False,
        )
    except Exception as e:
        for i in both_idx:
            rec = {"pdb": paths_chunk[i].name,
                   "error": f"batch_prepare: {type(e).__name__}: {e}"}
            out_fh.write(json.dumps(rec) + "\n")
            results.append(rec)
        out_fh.flush()
        return results

    for i, report in zip(both_idx, reports):
        rec = {"pdb": paths_chunk[i].name}
        try:
            if report.skipped_no_protein:
                rec["skipped"] = "no_protein"
            else:
                rec["initial_energy"] = float(report.initial_energy)
                rec["final_energy"] = float(report.final_energy)
                rec["minimizer_steps"] = int(report.minimizer_steps)
                rec["converged"] = bool(report.converged)
                ca_pre = ca_pre_by_idx[i]
                ca_post = proteon.extract_ca_coords(work_by_idx[i])
                rec["n_ca_pre"] = int(len(ca_pre))
                rec["n_ca_post"] = int(len(ca_post))
                if ca_pre.shape != ca_post.shape:
                    rec["error"] = f"CA shape mismatch {ca_pre.shape} vs {ca_post.shape}"
                else:
                    rec.update(tm_pair(ca_pre, ca_post))
        except Exception as e:
            rec["error"] = f"{type(e).__name__}: {e}"
        out_fh.write(json.dumps(rec) + "\n")
        results.append(rec)

    out_fh.flush()
    return results


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    pdb_list = os.environ.get("PROTEON_PDB_LIST")
    if pdb_list:
        list_path = Path(pdb_list)
        if not list_path.is_file():
            raise SystemExit(f"PDB list file not found: {list_path}")
        with open(list_path) as f:
            paths = [line.strip() for line in f if line.strip()]
        paths = [p for p in paths if Path(p).is_file()][:N]
        sample = [Path(p) for p in paths]
        print(f"Loaded {len(sample)} PDB paths from {list_path}", flush=True)
    else:
        pdbs = sorted(p.name for p in PDB_DIR.glob("*.pdb"))
        rng = random.Random(SEED)
        rng.shuffle(pdbs)
        sample = [PDB_DIR / name for name in pdbs[:N]]
        print(f"Sampled {len(sample)} PDBs from {PDB_DIR} (seed={SEED})", flush=True)

    done_names: set[str] = set()
    if OUT.exists():
        with open(OUT) as f:
            for line in f:
                try:
                    done_names.add(json.loads(line)["pdb"])
                except Exception:
                    pass
        print(f"Resuming: {len(done_names)} PDBs already in {OUT}", flush=True)
    pending = [p for p in sample if p.name not in done_names]
    print(f"Total sample: {len(sample)} PDBs; {len(pending)} pending after resume", flush=True)
    print(f"Writing to {OUT}", flush=True)

    t0 = time.perf_counter()
    n_ok = n_fail = n_skip = 0
    with open(OUT, "a") as f:
        for start in range(0, len(pending), CHUNK):
            chunk = pending[start:start + CHUNK]
            t_c = time.perf_counter()
            recs = process_chunk(chunk, f)
            dt = time.perf_counter() - t_c
            for r in recs:
                if "tm_score" in r:
                    n_ok += 1
                elif r.get("skipped"):
                    n_skip += 1
                else:
                    n_fail += 1
            done = start + len(chunk)
            progress = len(done_names) + done
            el = time.perf_counter() - t0
            rate = done / el if el > 0 else 0
            eta = (len(pending) - done) / rate if rate > 0 else 0
            print(
                f"[{progress}/{len(sample)}] chunk={len(chunk)} in {dt:.1f}s "
                f"({len(chunk)/dt:.2f}/s)  ok={n_ok} fail={n_fail} skip={n_skip}  "
                f"total_rate={rate:.2f}/s  eta={eta/60:.1f}min",
                flush=True,
            )

    elapsed = time.perf_counter() - t0
    print(f"\nDone. ok={n_ok} fail={n_fail} skip={n_skip} in {elapsed/60:.1f} min "
          f"({n_ok/elapsed:.2f} struct/s)", flush=True)

    tms, rmsds = [], []
    with open(OUT) as f:
        for l in f:
            r = json.loads(l)
            if "tm_score" in r:
                tms.append(r["tm_score"])
                rmsds.append(r["rmsd"])
    if tms:
        tms = np.array(tms); rmsds = np.array(rmsds)
        print(f"\nProteon AMBER96 TM (n={len(tms)}):")
        print(f"  mean={tms.mean():.4f}  median={np.median(tms):.4f}")
        print(f"  min={tms.min():.4f}  p05={np.percentile(tms,5):.4f}")
        print(f"RMSD: mean={rmsds.mean():.3f}  median={np.median(rmsds):.3f}  max={rmsds.max():.3f}")


if __name__ == "__main__":
    main()
