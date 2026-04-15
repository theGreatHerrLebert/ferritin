"""Validation and QA reporting for corpus releases."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from .corpus_release import CorpusReleaseManifest, load_corpus_release_manifest


@dataclass
class ValidationIssue:
    severity: str
    code: str
    message: str


@dataclass
class CorpusValidationReport:
    release_id: str
    artifact_type: str = "validation_report"
    format: str = "ferritin.corpus_validation.v0"
    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    ok: bool = True
    counts: Dict[str, int] = field(default_factory=dict)
    split_counts: Dict[str, int] = field(default_factory=dict)
    failure_breakdown: Dict[str, int] = field(default_factory=dict)
    completeness: Dict[str, float] = field(default_factory=dict)
    issues: List[ValidationIssue] = field(default_factory=list)


def validate_corpus_release(
    corpus_release_manifest: str | Path,
    *,
    out_path: str | Path | None = None,
) -> CorpusValidationReport:
    """Validate a corpus release and optionally write a JSON QA report."""
    manifest = load_corpus_release_manifest(corpus_release_manifest)
    report = CorpusValidationReport(
        release_id=manifest.release_id,
        counts={
            "prepared": manifest.count_prepared,
            "sequence_examples": manifest.count_sequence_examples,
            "structure_examples": manifest.count_structure_examples,
            "training_examples": manifest.count_training_examples,
        },
        split_counts=dict(manifest.split_counts),
        failure_breakdown=dict(manifest.failure_breakdown),
    )

    _check_count_consistency(manifest, report)
    _check_training_release(manifest, report)
    _check_structure_tensor_completeness(manifest, report)
    report.ok = not any(issue.severity == "error" for issue in report.issues)

    if out_path is not None:
        Path(out_path).write_text(
            json.dumps(
                {
                    **asdict(report),
                    "issues": [asdict(issue) for issue in report.issues],
                },
                indent=2,
            ),
            encoding="utf-8",
        )
    return report


def _check_count_consistency(manifest: CorpusReleaseManifest, report: CorpusValidationReport) -> None:
    if manifest.count_training_examples > manifest.count_sequence_examples:
        report.issues.append(
            ValidationIssue("error", "training_exceeds_sequence", "training example count exceeds sequence example count")
        )
    if manifest.count_training_examples > manifest.count_structure_examples:
        report.issues.append(
            ValidationIssue("error", "training_exceeds_structure", "training example count exceeds structure example count")
        )
    if manifest.count_prepared and manifest.count_structure_examples > manifest.count_prepared:
        report.issues.append(
            ValidationIssue("warning", "structure_exceeds_prepared", "structure example count exceeds prepared record count")
        )


def _check_training_release(manifest: CorpusReleaseManifest, report: CorpusValidationReport) -> None:
    if manifest.training_release is None:
        report.issues.append(ValidationIssue("warning", "missing_training_release", "no training release linked"))
        return
    rows = _load_jsonl(Path(manifest.training_release) / "training_examples.jsonl")
    if len(rows) != manifest.count_training_examples:
        report.issues.append(
            ValidationIssue("error", "training_row_count_mismatch", "training_examples.jsonl row count does not match manifest count")
        )
    split_counts: Dict[str, int] = {}
    duplicate_ids = set()
    seen = set()
    for row in rows:
        rid = row["record_id"]
        if rid in seen:
            duplicate_ids.add(rid)
        seen.add(rid)
        split = str(row.get("split", "train"))
        split_counts[split] = split_counts.get(split, 0) + 1
    if duplicate_ids:
        report.issues.append(
            ValidationIssue("error", "duplicate_training_ids", f"duplicate training record_ids found: {sorted(duplicate_ids)[:5]}")
        )
    if split_counts != manifest.split_counts:
        report.issues.append(
            ValidationIssue("error", "split_count_mismatch", "training split counts do not match corpus manifest")
        )


def _check_structure_tensor_completeness(manifest: CorpusReleaseManifest, report: CorpusValidationReport) -> None:
    if manifest.structure_release is None:
        report.issues.append(ValidationIssue("warning", "missing_structure_release", "no structure release linked"))
        return
    if manifest.count_structure_examples == 0:
        report.issues.append(
            ValidationIssue("warning", "no_structure_examples", "structure release contains no supervision examples")
        )
        return
    tensor_path = Path(manifest.structure_release) / "examples" / "tensors.npz"
    if not tensor_path.exists():
        report.issues.append(ValidationIssue("error", "missing_structure_tensors", "structure tensors.npz is missing"))
        return
    payload = np.load(tensor_path, allow_pickle=False)
    seq_mask = np.asarray(payload["seq_mask"], dtype=np.float32)
    rigid_exists = np.asarray(payload["rigidgroups_gt_exists"], dtype=np.float32)
    pseudo_mask = np.asarray(payload["pseudo_beta_mask"], dtype=np.float32)
    chi_mask = np.asarray(payload["chi_mask"], dtype=np.float32)

    valid_residues = float(seq_mask.sum())
    if valid_residues <= 0:
        report.issues.append(ValidationIssue("error", "no_valid_structure_residues", "structure seq_mask has no valid residues"))
        return
    report.completeness["pseudo_beta_fraction"] = float(pseudo_mask.sum() / valid_residues)
    report.completeness["rigidgroup_frame_fraction"] = float(rigid_exists.sum() / max(valid_residues * 8.0, 1.0))
    report.completeness["chi_angle_fraction"] = float(chi_mask.sum() / max(valid_residues * 4.0, 1.0))

    if report.completeness["pseudo_beta_fraction"] < 0.95:
        report.issues.append(
            ValidationIssue("warning", "low_pseudo_beta_fraction", "pseudo-beta completeness is below 95%")
        )


def _load_jsonl(path: Path) -> List[dict]:
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
