#!/usr/bin/env python3
"""Diagnose Foldseek-vs-ferritin retrieval benchmark outputs."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from statistics import mean
from typing import Any


def _source(row: dict[str, Any] | None) -> str | None:
    if not row:
        return None
    value = row.get("source_path")
    return None if value is None else str(value)


def _thresholds(data: dict[str, Any]) -> list[str]:
    configured = data.get("corpus", {}).get("thresholds")
    if configured:
        return [str(threshold) for threshold in configured]
    seen: set[str] = set()
    for row in data.get("per_query", []):
        seen.update(str(key) for key in row.get("thresholds", {}))
    return sorted(seen, key=float)


def _recall(row: dict[str, Any], threshold: str, method: str) -> float | None:
    threshold_row = row.get("thresholds", {}).get(threshold, {})
    key = "recall_at_k" if method == "foldseek" else "ferritin_recall_at_k"
    value = threshold_row.get(key)
    return None if value is None else float(value)


def _classify_top1(row: dict[str, Any]) -> str:
    truth = _source(row.get("best_truth_nonself"))
    foldseek = _source(row.get("foldseek_top_nonself"))
    ferritin = _source(row.get("ferritin_top_nonself"))
    if truth is None:
        return "no_truth"
    foldseek_exact = foldseek == truth
    ferritin_exact = ferritin == truth
    if foldseek_exact and ferritin_exact:
        return "both_exact"
    if foldseek_exact:
        return "foldseek_only"
    if ferritin_exact:
        return "ferritin_only"
    return "both_miss"


def _classify_recall(row: dict[str, Any], threshold: str, *, epsilon: float) -> str:
    foldseek = _recall(row, threshold, "foldseek")
    ferritin = _recall(row, threshold, "ferritin")
    if foldseek is None or ferritin is None:
        return "missing"
    delta = ferritin - foldseek
    if delta > epsilon:
        return "ferritin_better"
    if delta < -epsilon:
        return "foldseek_better"
    return "tie"


def _rank_of_source(hits: list[dict[str, Any]], source: str | None) -> int | None:
    if source is None:
        return None
    source_name = Path(source).name
    for idx, hit in enumerate(hits, start=1):
        if Path(str(hit.get("source_path", ""))).name == source_name:
            return idx
    return None


def summarize(data: dict[str, Any], *, primary_threshold: str | None = None, epsilon: float = 1e-9) -> dict[str, Any]:
    rows = list(data.get("per_query", []))
    thresholds = _thresholds(data)
    if not primary_threshold:
        primary_threshold = "0.7" if "0.7" in thresholds else (thresholds[len(thresholds) // 2] if thresholds else "")
    if primary_threshold and primary_threshold not in thresholds:
        raise ValueError(f"Threshold {primary_threshold!r} not present in benchmark output")

    top1_counts = Counter(_classify_top1(row) for row in rows)
    recall_counts = Counter(_classify_recall(row, primary_threshold, epsilon=epsilon) for row in rows) if primary_threshold else Counter()

    threshold_summary = {}
    for threshold in thresholds:
        foldseek_values = [_recall(row, threshold, "foldseek") for row in rows]
        ferritin_values = [_recall(row, threshold, "ferritin") for row in rows]
        paired = [
            (float(foldseek), float(ferritin))
            for foldseek, ferritin in zip(foldseek_values, ferritin_values)
            if foldseek is not None and ferritin is not None
        ]
        threshold_summary[threshold] = {
            "foldseek_mean": round(mean(foldseek for foldseek, _ in paired), 4) if paired else None,
            "ferritin_mean": round(mean(ferritin for _, ferritin in paired), 4) if paired else None,
            "ferritin_minus_foldseek": round(mean(ferritin - foldseek for foldseek, ferritin in paired), 4) if paired else None,
            "counts": dict(Counter(
                _classify_recall(row, threshold, epsilon=epsilon)
                for row in rows
            )),
        }

    query_diagnostics = []
    for row in rows:
        foldseek_recall = _recall(row, primary_threshold, "foldseek") if primary_threshold else None
        ferritin_recall = _recall(row, primary_threshold, "ferritin") if primary_threshold else None
        delta = None if foldseek_recall is None or ferritin_recall is None else ferritin_recall - foldseek_recall
        truth = _source(row.get("best_truth_nonself"))
        query_diagnostics.append({
            "query": row.get("query"),
            "truth": truth,
            "truth_tm_score": None if not row.get("best_truth_nonself") else row["best_truth_nonself"].get("tm_score"),
            "foldseek_top": _source(row.get("foldseek_top_nonself")),
            "ferritin_top": _source(row.get("ferritin_top_nonself")),
            "top1_class": _classify_top1(row),
            "foldseek_recall": foldseek_recall,
            "ferritin_recall": ferritin_recall,
            "ferritin_minus_foldseek": delta,
            "truth_rank_in_foldseek_trace": _rank_of_source(row.get("foldseek_top_hits", []), truth),
            "truth_rank_in_ferritin_trace": _rank_of_source(row.get("ferritin_top_hits", []), truth),
        })

    query_diagnostics.sort(
        key=lambda row: (
            0 if row["ferritin_minus_foldseek"] is not None else 1,
            row["ferritin_minus_foldseek"] if row["ferritin_minus_foldseek"] is not None else 0.0,
            str(row["query"]),
        )
    )

    skipped_truth = data.get("skipped_truth_candidates", {})
    truth_skip_metadata_known = bool(skipped_truth) or not data.get("corpus", {}).get("truth_cache_hit")
    return {
        "corpus": data.get("corpus", {}),
        "timing": data.get("timing", {}),
        "metrics": data.get("metrics", {}),
        "primary_threshold": primary_threshold,
        "n_queries": len(rows),
        "top1_counts": dict(top1_counts),
        "recall_counts": dict(recall_counts),
        "threshold_summary": threshold_summary,
        "worst_ferritin_queries": query_diagnostics[:10],
        "best_ferritin_queries": list(reversed(query_diagnostics[-10:])),
        "truth_skip_metadata_known": truth_skip_metadata_known,
        "n_skipped_truth_candidates": (
            sum(len(items) for items in skipped_truth.values())
            if truth_skip_metadata_known
            else None
        ),
        "n_queries_with_skipped_truth_candidates": (
            sum(1 for items in skipped_truth.values() if items)
            if truth_skip_metadata_known
            else None
        ),
        "skipped_ferritin_queries": data.get("skipped_ferritin_queries", {}),
        "has_candidate_traces": any(
            "foldseek_top_hits" in row or "ferritin_top_hits" in row
            for row in rows
        ),
    }


def render_markdown(summary: dict[str, Any]) -> str:
    def fmt(value: Any) -> str:
        if value is None:
            return ""
        if isinstance(value, float):
            return str(round(value, 4))
        return str(value)

    lines = [
        "# Foldseek vs ferritin Retrieval Diagnostics",
        "",
        "## Headline",
        "",
        f"- Queries: {summary['n_queries']}",
        f"- Primary threshold: TM >= {summary['primary_threshold']}",
        f"- Candidate traces present: {summary['has_candidate_traces']}",
        "- Skipped truth candidates: "
        + (
            f"{summary['n_skipped_truth_candidates']} across {summary['n_queries_with_skipped_truth_candidates']} queries"
            if summary["truth_skip_metadata_known"]
            else "unknown; benchmark reused a truth cache without skip metadata"
        ),
        f"- Skipped ferritin queries: {len(summary['skipped_ferritin_queries'])}",
        "",
        "## Recall By Threshold",
        "",
        "| TM threshold | Foldseek mean | ferritin mean | ferritin - Foldseek |",
        "|---:|---:|---:|---:|",
    ]
    for threshold, row in summary["threshold_summary"].items():
        lines.append(
            f"| {threshold} | {row['foldseek_mean']} | {row['ferritin_mean']} | {row['ferritin_minus_foldseek']} |"
        )

    lines.extend([
        "",
        "## Top-1 Classes",
        "",
    ])
    for key in ["both_exact", "foldseek_only", "ferritin_only", "both_miss", "no_truth"]:
        lines.append(f"- {key}: {summary['top1_counts'].get(key, 0)}")

    lines.extend([
        "",
        "## Primary Recall Classes",
        "",
    ])
    for key in ["foldseek_better", "ferritin_better", "tie", "missing"]:
        lines.append(f"- {key}: {summary['recall_counts'].get(key, 0)}")

    lines.extend([
        "",
        "## Worst ferritin Deltas",
        "",
        "| Query | Truth | Truth TM | Foldseek top | ferritin top | Foldseek recall | ferritin recall | Delta | ferritin truth rank |",
        "|---|---|---:|---|---|---:|---:|---:|---:|",
    ])
    for row in summary["worst_ferritin_queries"]:
        lines.append(
            "| {query} | {truth} | {truth_tm_score} | {foldseek_top} | {ferritin_top} | "
            "{foldseek_recall} | {ferritin_recall} | {ferritin_minus_foldseek} | {truth_rank_in_ferritin_trace} |".format(
                **{key: fmt(value) for key, value in row.items()}
            )
        )

    lines.extend([
        "",
        "## Best ferritin Deltas",
        "",
        "| Query | Truth | Truth TM | Foldseek top | ferritin top | Foldseek recall | ferritin recall | Delta | ferritin truth rank |",
        "|---|---|---:|---|---|---:|---:|---:|---:|",
    ])
    for row in summary["best_ferritin_queries"]:
        lines.append(
            "| {query} | {truth} | {truth_tm_score} | {foldseek_top} | {ferritin_top} | "
            "{foldseek_recall} | {ferritin_recall} | {ferritin_minus_foldseek} | {truth_rank_in_ferritin_trace} |".format(
                **{key: fmt(value) for key, value in row.items()}
            )
        )

    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Benchmark JSON from validation/bench_foldseek_retrieval.py")
    parser.add_argument("--threshold", default=None, help="Primary TM threshold for query-level deltas")
    parser.add_argument("--json-output", type=Path, default=None, help="Optional machine-readable diagnostics path")
    parser.add_argument("--markdown-output", type=Path, default=None, help="Optional markdown report path")
    args = parser.parse_args()

    data = json.loads(args.input.read_text(encoding="utf-8"))
    summary = summarize(data, primary_threshold=args.threshold)
    markdown = render_markdown(summary)
    print(markdown, end="")
    if args.json_output:
        args.json_output.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if args.markdown_output:
        args.markdown_output.write_text(markdown, encoding="utf-8")


if __name__ == "__main__":
    main()
