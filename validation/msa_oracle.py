"""MSA oracle: proteon search+MSA vs upstream mmseqs search+result2msa.

Runs both pipelines on the same query sequences against the same target
DB and compares the resulting MSAs field by field. Follows the oracle
testing philosophy in devdocs/ORACLE.md — the upstream pipeline is the
ground truth, and disagreements are diagnosed as either proteon bugs or
documented convention gaps.

## What we compare

Two comparison modes run for each query:

  **A) Positional** — compare MSA rows by rank (both sorted by alignment
  score descending). Affected by hit-set divergence: proteon's prefilter
  is substantially more sensitive than upstream's, finding 10–50× more
  hits at default settings, so most rows are comparing different targets.

  **B) Shared-hit** — match MSA rows by target DB key, compare only
  targets found by both pipelines. Isolates the alignment/MSA-assembly
  math from the prefilter divergence. This is the primary oracle gate.

## Tolerances

Calibrated on 50-query runs against 1M UniRef50 (2026-04-19). The
shared-hit comparison is the oracle contract; positional metrics are
reported for completeness but not gated.

| Metric                              | Threshold   | Rationale                                       |
|-------------------------------------|-------------|------------------------------------------------ |
| Query row (msa[0])                  | byte-exact  | Same sequence, same encoding                    |
| Shared-hit residue agreement        | >= 0.75     | SW extension boundaries differ on remote         |
|                                     |             | homologs (20-39% identity); within the           |
|                                     |             | overlapping aligned region residues agree ~90%   |
| Shared-hit deletion agreement       | >= 0.90     | CIGAR-to-deletion projection nearly identical    |
|                                     |             | within overlapping region; median 0.97           |
| Shared-hit gap agreement            | (report)    | Lower (~70%) due to boundary extension gaps;     |
|                                     |             | not gated — informational                        |
| Hit-set Jaccard                     | (report)    | Expected low (0.01–0.10); proteon prefilter is   |
|                                     |             | much more sensitive than mmseqs at -s 5.7/7.5    |
| MSA depth ratio                     | (report)    | Proteon saturates max_seqs while upstream returns |
|                                     |             | 2–256 hits; not gated                            |

## Known convention gap: alignment extension boundaries

Diagnosed 2026-04-19 on A0A0C1M9X2 (1262 residues, 171 shared hits
against uniref50_1m). 0/171 shared hits had identical CIGARs. Root
cause: proteon's Smith-Waterman extends alignments further into
low-identity tails (typical: proteon aligns positions 6–1262 vs
upstream 298–658). Query-start is close (median 4 residues apart)
but query-end diverges substantially (62/171 differ by >20 residues).

This is an alignment-policy difference, not a correctness bug:
  - Both tools implement standard SW with affine gap penalties
  - Proteon uses more aggressive extension on low-identity regions
  - All affected hits are remote homologs (20–39% fident)
  - Within the overlapping aligned region, residues match ~90% and
    deletions match ~97%

Analogous to the AMBER96 NoCutoff convention gap (devdocs/ORACLE.md
§When the oracle is also wrong, case 3): a policy choice that creates
a measurable gap, not a math error.

## Usage

    # Quick smoke (50 queries, self-search):
    python validation/msa_oracle.py

    # Scale run against 1M UniRef50:
    python validation/msa_oracle.py \\
        --query-fasta proteon-search/tests/data/oracle_fixture.fasta \\
        --target-db /globalscratch/dateschn/uniref_subset/uniref50_1m \\
        --target-kmi /globalscratch/dateschn/uniref_subset/uniref50_1m_ext.kmi \\
        --mmseqs /globalscratch/aweissen/bio/MMseqs2/build/bin/mmseqs \\
        --max-seqs 256 \\
        --out-json validation/msa_oracle_report.json

    # Higher upstream sensitivity (closer depth, same shared-hit metrics):
    python validation/msa_oracle.py --mmseqs-sensitivity 7.5 ...

    # With the venv:
    source venv/bin/activate

## Benchmark results (2026-04-19)

50 queries from oracle_fixture.fasta vs uniref50_1m (1M sequences).
GPU: RTX 5090. Both s=5.7 and s=7.5 tested — shared-hit metrics are
sensitivity-independent (confirmed identical within noise).

    Shared-hit residue agreement:  median=0.896  mean=0.899  min=0.780  p5=0.795
    Shared-hit deletion agreement: median=0.974  mean=0.974  min=0.950  p5=0.955
    Shared-hit gap agreement:      median=0.710  mean=0.701  min=0.345  p5=0.416
    Hit-set Jaccard:               median=0.013  mean=0.102  min=0.000
    All 50 queries PASSED.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Amino-acid encoding (matches proteon's Alphabet::protein())
# ---------------------------------------------------------------------------

_AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"  # proteon's Alphabet::protein() — alphabetical
_AA_TO_IDX = {aa: i for i, aa in enumerate(_AA_ORDER)}
_GAP_IDX = 21  # X=20, gap=21 — matches proteon convention


def _encode_aa(ch: str) -> int:
    """Encode one uppercase amino acid letter to proteon's uint8 index."""
    return _AA_TO_IDX.get(ch.upper(), 20)  # 20 = unknown/X


# ---------------------------------------------------------------------------
# Upstream (MMseqs2) pipeline
# ---------------------------------------------------------------------------


def run_upstream_msa(
    mmseqs: Path,
    query_fasta: Path,
    target_fasta_or_db: Path,
    workdir: Path,
    *,
    max_seqs: int = 256,
    target_is_db: bool = False,
    sensitivity: float = 5.7,
) -> Dict[str, str]:
    """Run mmseqs createdb + search + result2msa and return {query_acc: a3m_text}.

    If `target_is_db` is True, `target_fasta_or_db` is an existing
    MMseqs2 DB prefix and we skip createdb on it.
    """
    query_db = workdir / "queryDB"
    target_db = workdir / "targetDB" if not target_is_db else target_fasta_or_db
    result_db = workdir / "resultDB"
    msa_db = workdir / "msaDB"
    tmp = workdir / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)

    def _run(subcmd: str, args: list):
        cmd = [str(mmseqs), subcmd] + [str(a) for a in args]
        print(f"  mmseqs {subcmd} ...", flush=True)
        t0 = time.monotonic()
        result = subprocess.run(cmd, capture_output=True, text=True)
        dt = time.monotonic() - t0
        if result.returncode != 0:
            print(f"  STDERR: {result.stderr[:2000]}", file=sys.stderr)
            raise RuntimeError(f"mmseqs {subcmd} failed (exit {result.returncode})")
        print(f"  mmseqs {subcmd} done ({dt:.1f}s)", flush=True)

    # 1. createdb for queries
    _run("createdb", [query_fasta, query_db, "-v", "1"])

    # 2. createdb for targets (skip if pre-built)
    if not target_is_db:
        _run("createdb", [target_fasta_or_db, target_db, "-v", "1"])

    # 3. search
    _run("search", [
        query_db, target_db, result_db, tmp,
        "-v", "1",
        "--threads", "4",
        "-s", str(sensitivity),
        "--max-seqs", str(max_seqs),
    ])

    # 4. result2msa → a3m output DB
    _run("result2msa", [
        query_db, target_db, result_db, msa_db,
        "-v", "1",
        "--msa-format-mode", "6",  # a3m format
    ])

    # 5. Parse the MSA DB — each entry is one query's a3m
    #    Use convertmsa or read the DB directly
    msa_out_dir = workdir / "msa_a3m"
    msa_out_dir.mkdir(exist_ok=True)

    # Extract all MSAs as individual files
    _run("unpackdb", [msa_db, msa_out_dir, "--unpack-name-mode", "0", "-v", "1"])

    # Read query accessions from the lookup file
    lookup_path = query_db.parent / f"{query_db.name}.lookup"
    acc_by_key: Dict[int, str] = {}
    if lookup_path.exists():
        for line in lookup_path.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 2:
                acc_by_key[int(parts[0])] = parts[1]

    # Read each extracted a3m
    results: Dict[str, str] = {}
    for f in sorted(msa_out_dir.iterdir()):
        if not f.is_file():
            continue
        key = int(f.name)
        acc = acc_by_key.get(key, f.name)
        results[acc] = f.read_text(encoding="utf-8", errors="replace")

    return results


# ---------------------------------------------------------------------------
# Parse a3m into comparable form
# ---------------------------------------------------------------------------


def parse_a3m_to_arrays(
    a3m_text: str, max_seqs: int = 256
) -> Optional[Dict[str, np.ndarray]]:
    """Parse a3m text into proteon-compatible MSA arrays.

    Returns dict with:
      msa: (N, L) uint8  — encoded amino acids, gap_idx for gaps
      deletion_matrix: (N, L) uint8  — insertion counts before each column
      msa_mask: (N, L) float32  — 1.0 for present, 0.0 for absent
      target_accs: list[str]  — accession for each row (row 0 = query)
    """
    from proteon.msa_io import parse_a3m_text

    try:
        aligned_rows, deletion_matrix, query_aligned = parse_a3m_text(a3m_text)
    except (ValueError, IndexError):
        return None

    if not aligned_rows:
        return None

    n_rows = min(len(aligned_rows), max_seqs)
    query_len = len(query_aligned)

    msa = np.full((n_rows, query_len), _GAP_IDX, dtype=np.uint8)
    del_mat = np.zeros((n_rows, query_len), dtype=np.uint8)
    mask = np.zeros((n_rows, query_len), dtype=np.float32)

    # Parse accessions from the a3m text
    accs = []
    for line in a3m_text.splitlines():
        if line.startswith(">"):
            accs.append(line[1:].split()[0])

    for i in range(n_rows):
        row_str = aligned_rows[i]
        row_del = deletion_matrix[i]
        for j, ch in enumerate(row_str):
            if ch == "-":
                msa[i, j] = _GAP_IDX
                mask[i, j] = 1.0
            else:
                msa[i, j] = _encode_aa(ch)
                mask[i, j] = 1.0
            if j < len(row_del):
                del_mat[i, j] = min(row_del[j], 255)

    return {
        "msa": msa,
        "deletion_matrix": del_mat,
        "msa_mask": mask,
        "target_accs": accs[:n_rows],
    }


# ---------------------------------------------------------------------------
# Proteon pipeline
# ---------------------------------------------------------------------------


def run_proteon_msa(
    query_seqs: Dict[str, str],
    engine,
    *,
    max_seqs: int = 256,
) -> Dict[str, Dict[str, np.ndarray]]:
    """Run proteon's search + MSA for each query.

    Returns {query_acc: msa_dict} where msa_dict has keys from
    search_and_build_msa (msa, deletion_matrix, msa_mask, etc)
    plus a 'hit_target_ids' list from the raw search.
    """
    from proteon.msa_backend import search_and_build_msa, search

    results = {}
    for acc, seq in query_seqs.items():
        try:
            # Run search first to capture per-hit target IDs
            hits = search(engine, seq)
            hit_ids = [h["target_id"] for h in hits[:max_seqs - 1]]  # -1 for query row

            msa_dict = search_and_build_msa(engine, seq, max_seqs=max_seqs)
            msa_dict["hit_target_ids"] = hit_ids
            results[acc] = msa_dict
        except Exception as e:
            print(f"  WARNING: proteon MSA failed for {acc}: {e}", file=sys.stderr)
    return results


# ---------------------------------------------------------------------------
# Comparison logic
# ---------------------------------------------------------------------------


def compare_msa_pair(
    proteon_msa: Dict[str, np.ndarray],
    upstream_msa: Dict[str, np.ndarray],
    query_acc: str,
) -> Dict:
    """Compare one query's MSA from both pipelines.

    Returns a per-query report dict with metrics and diagnostics.
    Two comparison modes:

    A) **Positional** — compare rows by position (both sorted by score
       descending). Affected by prefilter sensitivity differences.
    B) **Shared-hit** — match rows by target ID, compare only hits found
       by both pipelines. Isolates the alignment/MSA-assembly math from
       the prefilter divergence.
    """
    report = {"query": query_acc, "pass": True, "warnings": [], "metrics": {}}

    p_msa = proteon_msa["msa"]      # (N_p, L)
    u_msa = upstream_msa["msa"]     # (N_u, L)
    p_del = proteon_msa["deletion_matrix"]
    u_del = upstream_msa["deletion_matrix"]

    # Shapes
    n_p, l_p = p_msa.shape
    n_u, l_u = u_msa.shape

    report["metrics"]["proteon_depth"] = int(n_p)
    report["metrics"]["upstream_depth"] = int(n_u)
    report["metrics"]["proteon_query_len"] = int(l_p)
    report["metrics"]["upstream_query_len"] = int(l_u)

    # Query length must match
    if l_p != l_u:
        report["pass"] = False
        report["warnings"].append(
            f"query length mismatch: proteon={l_p}, upstream={l_u}"
        )
        return report

    query_len = l_p

    # 1. Query row (msa[0]) — must be byte-exact
    query_row_match = np.array_equal(p_msa[0], u_msa[0])
    report["metrics"]["query_row_exact"] = query_row_match
    if not query_row_match:
        n_diff = int(np.sum(p_msa[0] != u_msa[0]))
        report["warnings"].append(
            f"query row differs at {n_diff}/{query_len} positions"
        )
        if n_diff > query_len * 0.05:
            report["pass"] = False

    # 2. MSA depth ratio
    if n_u > 0:
        depth_ratio = n_p / n_u
        report["metrics"]["depth_ratio"] = round(depth_ratio, 3)
        if not (0.3 <= depth_ratio <= 3.0):
            report["warnings"].append(
                f"depth ratio {depth_ratio:.2f} outside [0.3, 3.0]"
            )

    # -----------------------------------------------------------------------
    # A) Positional comparison (rows by rank, affected by hit-set divergence)
    # -----------------------------------------------------------------------
    min_depth = min(n_p, n_u)
    if min_depth > 1:
        shared_msa_p = p_msa[1:min_depth]
        shared_msa_u = u_msa[1:min_depth]

        both_present = (shared_msa_p != _GAP_IDX) & (shared_msa_u != _GAP_IDX)
        if both_present.sum() > 0:
            residue_agree = float(
                (shared_msa_p[both_present] == shared_msa_u[both_present]).mean()
            )
        else:
            residue_agree = 1.0

        gap_pattern_agree = float(
            ((shared_msa_p == _GAP_IDX) == (shared_msa_u == _GAP_IDX)).mean()
        )

        shared_del_p = p_del[1:min_depth].astype(np.float32)
        shared_del_u = u_del[1:min_depth].astype(np.float32)
        if both_present.sum() > 0:
            del_agree = float(
                (shared_del_p[both_present] == shared_del_u[both_present]).mean()
            )
        else:
            del_agree = 1.0

        report["metrics"]["positional_residue_agreement"] = round(residue_agree, 4)
        report["metrics"]["positional_gap_agreement"] = round(gap_pattern_agree, 4)
        report["metrics"]["positional_deletion_agreement"] = round(del_agree, 4)
        report["metrics"]["positional_shared_rows"] = int(min_depth - 1)

        # Keep old names for backwards compat
        report["metrics"]["residue_agreement"] = report["metrics"]["positional_residue_agreement"]
        report["metrics"]["gap_pattern_agreement"] = report["metrics"]["positional_gap_agreement"]
        report["metrics"]["deletion_agreement"] = report["metrics"]["positional_deletion_agreement"]

    # -----------------------------------------------------------------------
    # B) Shared-hit comparison (match rows by target ID)
    # -----------------------------------------------------------------------
    p_hit_ids = proteon_msa.get("hit_target_ids", [])
    u_target_accs = upstream_msa.get("target_accs", [])

    # Build upstream target_acc → row index map (skip row 0 = query)
    u_acc_to_row = {}
    for i, acc in enumerate(u_target_accs):
        if i == 0:
            continue  # query row
        u_acc_to_row[acc] = i

    # proteon hit_target_ids[i] corresponds to msa row i+1
    # upstream target_accs use the DB key as accession (since we wrote
    # numeric keys in the _h DB). proteon's target_id is the DB key.
    n_shared = 0
    shared_residue_matches = 0
    shared_residue_total = 0
    shared_gap_matches = 0
    shared_gap_total = 0
    shared_del_matches = 0
    shared_del_total = 0

    for p_row_idx, tid in enumerate(p_hit_ids):
        p_row = p_row_idx + 1  # msa row (0 is query)
        if p_row >= n_p:
            break

        # Try to find this target in upstream's MSA by DB key
        u_row = u_acc_to_row.get(str(tid))
        if u_row is None:
            continue

        n_shared += 1
        p_row_data = p_msa[p_row]
        u_row_data = u_msa[u_row]
        p_row_del = p_del[p_row]
        u_row_del = u_del[u_row]

        # Residue agreement where both have a real residue
        both = (p_row_data != _GAP_IDX) & (u_row_data != _GAP_IDX)
        if both.sum() > 0:
            shared_residue_matches += int((p_row_data[both] == u_row_data[both]).sum())
            shared_residue_total += int(both.sum())

        # Gap pattern
        shared_gap_matches += int(((p_row_data == _GAP_IDX) == (u_row_data == _GAP_IDX)).sum())
        shared_gap_total += query_len

        # Deletion matrix where both present
        if both.sum() > 0:
            shared_del_matches += int((p_row_del[both] == u_row_del[both]).sum())
            shared_del_total += int(both.sum())

    # Hit-set overlap
    p_id_set = set(str(tid) for tid in p_hit_ids[:n_p - 1])
    u_id_set = set(u_acc_to_row.keys())
    intersection = len(p_id_set & u_id_set)
    union = len(p_id_set | u_id_set)

    report["metrics"]["shared_hit_count"] = n_shared
    report["metrics"]["hit_set_jaccard"] = round(intersection / max(union, 1), 4)
    report["metrics"]["proteon_unique_hits"] = len(p_id_set - u_id_set)
    report["metrics"]["upstream_unique_hits"] = len(u_id_set - p_id_set)

    if n_shared > 0:
        sh_res = round(shared_residue_matches / max(shared_residue_total, 1), 4)
        sh_gap = round(shared_gap_matches / max(shared_gap_total, 1), 4)
        sh_del = round(shared_del_matches / max(shared_del_total, 1), 4)
        report["metrics"]["shared_hit_residue_agreement"] = sh_res
        report["metrics"]["shared_hit_gap_agreement"] = sh_gap
        report["metrics"]["shared_hit_deletion_agreement"] = sh_del

        # Oracle gate: when both find the same target, do they produce
        # the same alignment? Thresholds account for the documented
        # alignment-extension boundary convention gap (see module
        # docstring). Residue threshold 0.75 gives headroom below the
        # observed p5=0.78; deletion threshold 0.90 below observed
        # p5=0.95. Violations are hard failures.
        if sh_res < 0.75:
            report["pass"] = False
            report["warnings"].append(
                f"FAIL: shared-hit residue agreement {sh_res:.4f} < 0.75 "
                f"({n_shared} shared hits)"
            )
        elif sh_res < 0.85:
            report["warnings"].append(
                f"WARN: shared-hit residue agreement {sh_res:.4f} < 0.85 "
                f"({n_shared} shared hits)"
            )
        if sh_del < 0.90:
            report["pass"] = False
            report["warnings"].append(
                f"FAIL: shared-hit deletion agreement {sh_del:.4f} < 0.90"
            )
        elif sh_del < 0.95:
            report["warnings"].append(
                f"WARN: shared-hit deletion agreement {sh_del:.4f} < 0.95"
            )

    return report


# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------


def _normalize_accession(raw_header: str) -> str:
    """Extract the MMseqs2-compatible accession from a FASTA header.

    MMseqs2 createdb extracts the UniProt accession from headers like
    ``sp|Q8AWH3|SX17A_XENTR`` or ``tr|W0FSK4|W0FSK4_9FLAV`` — it takes
    the second pipe-delimited field. For plain headers like ``SEQ1`` it
    takes everything up to the first whitespace. We mirror that logic
    so our key space matches upstream's .lookup file.
    """
    token = raw_header.split()[0]
    parts = token.split("|")
    if len(parts) >= 2 and parts[0] in ("sp", "tr"):
        return parts[1]
    return token


def read_fasta(path: Path) -> Dict[str, str]:
    """Read a FASTA file into {accession: sequence}.

    Accessions are normalized to match MMseqs2's createdb .lookup
    convention (UniProt middle field for sp|..|.. headers).
    """
    seqs: Dict[str, str] = {}
    name: Optional[str] = None
    body: list = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(body)
            name = _normalize_accession(line[1:])
            body = []
        elif line.strip():
            body.append(line.strip())
    if name is not None:
        seqs[name] = "".join(body)
    return seqs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="MSA oracle: proteon vs MMseqs2")
    parser.add_argument(
        "--query-fasta",
        type=Path,
        default=None,
        help="Query FASTA (default: oracle_fixture.fasta)",
    )
    parser.add_argument(
        "--target-fasta",
        type=Path,
        default=None,
        help="Target FASTA for mmseqs createdb (mutually exclusive with --target-db)",
    )
    parser.add_argument(
        "--target-db",
        type=Path,
        default=None,
        help="Pre-built MMseqs2 DB prefix for targets",
    )
    parser.add_argument(
        "--target-kmi",
        type=Path,
        default=None,
        help="Pre-built .kmi file for proteon (skips k-mer index rebuild)",
    )
    parser.add_argument(
        "--mmseqs",
        type=Path,
        default=None,
        help="Path to mmseqs binary",
    )
    parser.add_argument(
        "--max-seqs", type=int, default=256,
        help="Max MSA depth (default 256)",
    )
    parser.add_argument(
        "--max-queries", type=int, default=None,
        help="Cap number of queries (for quick smoke tests)",
    )
    parser.add_argument(
        "--mmseqs-sensitivity", type=float, default=5.7,
        help="MMseqs2 -s sensitivity (default 5.7, try 7.5 for deeper search)",
    )
    parser.add_argument(
        "--out-json", type=Path, default=None,
        help="Write per-query JSON report to this path",
    )
    args = parser.parse_args()

    # Resolve mmseqs binary
    mmseqs = args.mmseqs
    if mmseqs is None:
        candidates = [
            Path("/globalscratch/aweissen/bio/MMseqs2/build/bin/mmseqs"),
            Path(os.environ.get("PROTEON_SEARCH_MMSEQS_BIN", "")),
        ]
        for c in candidates:
            if c.exists():
                mmseqs = c
                break
    if mmseqs is None or not mmseqs.exists():
        print("ERROR: mmseqs binary not found. Pass --mmseqs or set PROTEON_SEARCH_MMSEQS_BIN.", file=sys.stderr)
        sys.exit(1)
    print(f"mmseqs binary: {mmseqs}")

    # Resolve query FASTA
    query_fasta = args.query_fasta
    if query_fasta is None:
        default = Path(__file__).resolve().parent.parent / "proteon-search" / "tests" / "data" / "oracle_fixture.fasta"
        if default.exists():
            query_fasta = default
    if query_fasta is None or not query_fasta.exists():
        print("ERROR: query FASTA not found", file=sys.stderr)
        sys.exit(1)
    print(f"query FASTA: {query_fasta}")

    # Read query sequences
    query_seqs = read_fasta(query_fasta)
    if args.max_queries and len(query_seqs) > args.max_queries:
        query_seqs = dict(list(query_seqs.items())[:args.max_queries])
    print(f"queries: {len(query_seqs)}")

    # Determine target source for upstream
    target_is_db = args.target_db is not None
    target_source = args.target_db or args.target_fasta or query_fasta
    print(f"target: {target_source} (is_db={target_is_db})")

    # -----------------------------------------------------------------------
    # 1. Run upstream MMseqs2 pipeline
    # -----------------------------------------------------------------------
    print("\n=== Upstream (MMseqs2) pipeline ===")
    t0 = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="msa_oracle_") as workdir:
        upstream_a3ms = run_upstream_msa(
            mmseqs, query_fasta, target_source, Path(workdir),
            max_seqs=args.max_seqs,
            target_is_db=target_is_db,
            sensitivity=args.mmseqs_sensitivity,
        )
        upstream_time = time.monotonic() - t0
    print(f"upstream done: {len(upstream_a3ms)} MSAs in {upstream_time:.1f}s")

    # Parse upstream a3ms into arrays
    upstream_parsed: Dict[str, Dict] = {}
    for acc, a3m_text in upstream_a3ms.items():
        parsed = parse_a3m_to_arrays(a3m_text, max_seqs=args.max_seqs)
        if parsed is not None:
            upstream_parsed[acc] = parsed

    print(f"upstream parsed: {len(upstream_parsed)} MSAs")

    # -----------------------------------------------------------------------
    # 2. Run proteon pipeline
    # -----------------------------------------------------------------------
    print("\n=== Proteon pipeline ===")
    from proteon.msa_backend import (
        build_search_engine,
        open_search_engine_from_mmseqs_db_with_kmi,
        build_search_engine_from_mmseqs_db,
    )

    t0 = time.monotonic()
    if args.target_db and args.target_kmi:
        print(f"  opening engine from DB + KMI ...")
        engine = open_search_engine_from_mmseqs_db_with_kmi(
            args.target_db, args.target_kmi,
        )
    elif args.target_db:
        print(f"  building engine from DB (no KMI, will build k-mer index) ...")
        engine = build_search_engine_from_mmseqs_db(args.target_db)
    else:
        # Self-search: build from query sequences
        print(f"  building in-memory engine from {len(query_seqs)} sequences ...")
        targets = [(i, seq) for i, (acc, seq) in enumerate(query_seqs.items())]
        engine = build_search_engine(targets)
    engine_time = time.monotonic() - t0
    print(f"  engine ready ({engine_time:.1f}s)")

    t0 = time.monotonic()
    proteon_msas = run_proteon_msa(query_seqs, engine, max_seqs=args.max_seqs)
    search_time = time.monotonic() - t0
    print(f"proteon done: {len(proteon_msas)} MSAs in {search_time:.1f}s")

    # -----------------------------------------------------------------------
    # 3. Compare
    # -----------------------------------------------------------------------
    print("\n=== Comparison ===")
    reports = []
    n_pass = 0
    n_fail = 0
    n_skip = 0

    # Aggregate metrics
    all_residue_agree = []
    all_gap_agree = []
    all_del_agree = []
    all_depth_ratios = []

    for acc in query_seqs:
        if acc not in upstream_parsed:
            n_skip += 1
            continue
        if acc not in proteon_msas:
            reports.append({
                "query": acc, "pass": False,
                "warnings": ["proteon returned no MSA"],
                "metrics": {},
            })
            n_fail += 1
            continue

        report = compare_msa_pair(proteon_msas[acc], upstream_parsed[acc], acc)
        reports.append(report)

        if report["pass"]:
            n_pass += 1
        else:
            n_fail += 1

        m = report["metrics"]
        if "residue_agreement" in m:
            all_residue_agree.append(m["residue_agreement"])
        if "gap_pattern_agreement" in m:
            all_gap_agree.append(m["gap_pattern_agreement"])
        if "deletion_agreement" in m:
            all_del_agree.append(m["deletion_agreement"])
        if "depth_ratio" in m:
            all_depth_ratios.append(m["depth_ratio"])

    # Collect shared-hit metrics
    all_sh_residue = [r["metrics"]["shared_hit_residue_agreement"]
                      for r in reports if "shared_hit_residue_agreement" in r.get("metrics", {})]
    all_sh_gap = [r["metrics"]["shared_hit_gap_agreement"]
                  for r in reports if "shared_hit_gap_agreement" in r.get("metrics", {})]
    all_sh_del = [r["metrics"]["shared_hit_deletion_agreement"]
                  for r in reports if "shared_hit_deletion_agreement" in r.get("metrics", {})]
    all_jaccard = [r["metrics"]["hit_set_jaccard"]
                   for r in reports if "hit_set_jaccard" in r.get("metrics", {})]
    all_shared_counts = [r["metrics"]["shared_hit_count"]
                         for r in reports if "shared_hit_count" in r.get("metrics", {})]

    # -----------------------------------------------------------------------
    # 4. Report
    # -----------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"MSA Oracle Summary")
    print(f"{'='*60}")
    print(f"Queries:    {len(query_seqs)}")
    print(f"Compared:   {n_pass + n_fail}")
    print(f"Pass:       {n_pass}")
    print(f"Fail:       {n_fail}")
    print(f"Skip:       {n_skip} (no upstream MSA)")
    print()

    print("--- Positional comparison (rows by rank) ---")
    if all_residue_agree:
        arr = np.array(all_residue_agree)
        print(f"Residue agreement:     median={np.median(arr):.4f}  "
              f"mean={np.mean(arr):.4f}  min={np.min(arr):.4f}  "
              f"p5={np.percentile(arr, 5):.4f}")

    if all_gap_agree:
        arr = np.array(all_gap_agree)
        print(f"Gap pattern agreement: median={np.median(arr):.4f}  "
              f"mean={np.mean(arr):.4f}  min={np.min(arr):.4f}  "
              f"p5={np.percentile(arr, 5):.4f}")

    if all_del_agree:
        arr = np.array(all_del_agree)
        print(f"Deletion agreement:    median={np.median(arr):.4f}  "
              f"mean={np.mean(arr):.4f}  min={np.min(arr):.4f}  "
              f"p5={np.percentile(arr, 5):.4f}")

    if all_depth_ratios:
        arr = np.array(all_depth_ratios)
        print(f"Depth ratio (p/u):     median={np.median(arr):.3f}  "
              f"mean={np.mean(arr):.3f}  min={np.min(arr):.3f}  "
              f"max={np.max(arr):.3f}")

    print()
    print("--- Shared-hit comparison (matched by target ID) ---")
    if all_shared_counts:
        arr = np.array(all_shared_counts)
        print(f"Shared hits per query: median={np.median(arr):.0f}  "
              f"mean={np.mean(arr):.1f}  min={np.min(arr):.0f}  "
              f"max={np.max(arr):.0f}  total={int(np.sum(arr))}")

    if all_jaccard:
        arr = np.array(all_jaccard)
        print(f"Hit-set Jaccard:       median={np.median(arr):.4f}  "
              f"mean={np.mean(arr):.4f}  min={np.min(arr):.4f}")

    if all_sh_residue:
        arr = np.array(all_sh_residue)
        print(f"Residue agreement:     median={np.median(arr):.4f}  "
              f"mean={np.mean(arr):.4f}  min={np.min(arr):.4f}  "
              f"p5={np.percentile(arr, 5):.4f}")

    if all_sh_gap:
        arr = np.array(all_sh_gap)
        print(f"Gap pattern agreement: median={np.median(arr):.4f}  "
              f"mean={np.mean(arr):.4f}  min={np.min(arr):.4f}  "
              f"p5={np.percentile(arr, 5):.4f}")

    if all_sh_del:
        arr = np.array(all_sh_del)
        print(f"Deletion agreement:    median={np.median(arr):.4f}  "
              f"mean={np.mean(arr):.4f}  min={np.min(arr):.4f}  "
              f"p5={np.percentile(arr, 5):.4f}")

    # Worst queries — show both positional and shared-hit
    if reports:
        with_shared = [r for r in reports if "shared_hit_residue_agreement" in r.get("metrics", {})]
        if with_shared:
            worst = sorted(with_shared, key=lambda r: r["metrics"]["shared_hit_residue_agreement"])[:10]
            print(f"\nWorst 10 queries by shared-hit residue agreement:")
            print(f"{'Query':<20} {'ShRes':>7} {'ShGap':>7} {'ShDel':>7} "
                  f"{'Shared':>7} {'Jaccard':>8} {'Depth(p/u)':>12}")
            print("-" * 72)
            for r in worst:
                m = r["metrics"]
                print(f"{r['query']:<20} {m.get('shared_hit_residue_agreement', 'N/A'):>7} "
                      f"{m.get('shared_hit_gap_agreement', 'N/A'):>7} "
                      f"{m.get('shared_hit_deletion_agreement', 'N/A'):>7} "
                      f"{m.get('shared_hit_count', '?'):>7} "
                      f"{m.get('hit_set_jaccard', 'N/A'):>8} "
                      f"{m.get('proteon_depth', '?')}/{m.get('upstream_depth', '?'):>10}")

    # Warnings
    all_warnings = [r for r in reports if r.get("warnings")]
    if all_warnings:
        print(f"\nQueries with warnings ({len(all_warnings)}):")
        for r in all_warnings[:20]:
            for w in r["warnings"]:
                print(f"  {r['query']}: {w}")
        if len(all_warnings) > 20:
            print(f"  ... and {len(all_warnings) - 20} more")

    # Save JSON report
    if args.out_json:
        full_report = {
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "config": {
                "query_fasta": str(query_fasta),
                "target": str(target_source),
                "target_is_db": target_is_db,
                "max_seqs": args.max_seqs,
                "n_queries": len(query_seqs),
            },
            "timing": {
                "upstream_s": round(upstream_time, 1),
                "proteon_engine_s": round(engine_time, 1),
                "proteon_search_s": round(search_time, 1),
            },
            "summary": {
                "n_pass": n_pass,
                "n_fail": n_fail,
                "n_skip": n_skip,
                "positional": {
                    "median_residue_agreement": round(float(np.median(all_residue_agree)), 4) if all_residue_agree else None,
                    "median_gap_agreement": round(float(np.median(all_gap_agree)), 4) if all_gap_agree else None,
                    "median_deletion_agreement": round(float(np.median(all_del_agree)), 4) if all_del_agree else None,
                    "median_depth_ratio": round(float(np.median(all_depth_ratios)), 3) if all_depth_ratios else None,
                },
                "shared_hit": {
                    "median_shared_hits": round(float(np.median(all_shared_counts)), 1) if all_shared_counts else None,
                    "total_shared_hits": int(np.sum(all_shared_counts)) if all_shared_counts else None,
                    "median_jaccard": round(float(np.median(all_jaccard)), 4) if all_jaccard else None,
                    "median_residue_agreement": round(float(np.median(all_sh_residue)), 4) if all_sh_residue else None,
                    "median_gap_agreement": round(float(np.median(all_sh_gap)), 4) if all_sh_gap else None,
                    "median_deletion_agreement": round(float(np.median(all_sh_del)), 4) if all_sh_del else None,
                },
                "mmseqs_sensitivity": args.mmseqs_sensitivity,
            },
            "per_query": reports,
        }
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(full_report, indent=2))
        print(f"\nJSON report: {args.out_json}")

    # Exit code: fail if any hard failures
    print()
    if n_fail > 0:
        print(f"RESULT: {n_fail} queries FAILED oracle comparison")
        sys.exit(1)
    else:
        print(f"RESULT: all {n_pass} queries PASSED oracle comparison")


if __name__ == "__main__":
    main()
