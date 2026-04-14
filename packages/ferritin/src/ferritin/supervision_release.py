"""Release-oriented wrappers around structure supervision exports."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from .supervision import StructureSupervisionExample
from .supervision_export import (
    SUPERVISION_EXPORT_FORMAT,
    export_structure_supervision_examples,
)


@dataclass
class FailureRecord:
    """Structured failure row for supervision/release pipelines."""

    record_id: str
    artifact_type: str = "failure_record"
    stage: str = "structure_supervision"
    status: str = "failed"
    failure_class: str = "internal_pipeline_error"
    message: str = ""
    source_id: Optional[str] = None
    prep_run_id: Optional[str] = None
    code_rev: Optional[str] = None
    config_rev: Optional[str] = None
    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    provenance: Dict[str, object] = field(default_factory=dict)


@dataclass
class StructureSupervisionReleaseManifest:
    """Machine-readable manifest for one supervision dataset release."""

    release_id: str
    artifact_type: str = "release_manifest"
    format: str = SUPERVISION_EXPORT_FORMAT
    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    code_rev: Optional[str] = None
    config_rev: Optional[str] = None
    count_examples: int = 0
    count_failures: int = 0
    example_export_dir: str = "examples"
    examples_file: str = "examples/examples.jsonl"
    tensor_file: str = "examples/tensors.npz"
    failure_file: str = "failures.jsonl"
    lengths: Dict[str, float] = field(default_factory=dict)
    sequence_lengths: List[int] = field(default_factory=list)
    provenance: Dict[str, object] = field(default_factory=dict)


def build_structure_supervision_release(
    examples: Iterable[StructureSupervisionExample],
    out_dir: str | Path,
    *,
    release_id: str,
    failures: Optional[Iterable[FailureRecord]] = None,
    code_rev: Optional[str] = None,
    config_rev: Optional[str] = None,
    provenance: Optional[Dict[str, object]] = None,
    overwrite: bool = False,
) -> Path:
    """Write a supervision release directory with examples, failures, and manifest."""
    example_list = list(examples)
    failure_list = list(failures or [])

    root = Path(out_dir)
    if root.exists() and not overwrite:
        raise FileExistsError(f"{root} already exists")
    root.mkdir(parents=True, exist_ok=True)

    example_dir = root / "examples"
    export_structure_supervision_examples(example_list, example_dir, overwrite=True)

    failure_path = root / "failures.jsonl"
    with failure_path.open("w", encoding="utf-8") as handle:
        for failure in failure_list:
            handle.write(json.dumps(asdict(failure), separators=(",", ":")))
            handle.write("\n")

    lengths = [ex.length for ex in example_list]
    manifest = StructureSupervisionReleaseManifest(
        release_id=release_id,
        code_rev=code_rev,
        config_rev=config_rev,
        count_examples=len(example_list),
        count_failures=len(failure_list),
        lengths=_length_summary(lengths),
        sequence_lengths=lengths,
        provenance=dict(provenance or {}),
    )
    (root / "release_manifest.json").write_text(
        json.dumps(asdict(manifest), indent=2),
        encoding="utf-8",
    )
    return root


def load_failure_records(path: str | Path) -> List[FailureRecord]:
    """Load failure records from a JSONL file."""
    rows = [
        json.loads(line)
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    return [FailureRecord(**row) for row in rows]


def _length_summary(lengths: List[int]) -> Dict[str, float]:
    if not lengths:
        return {"min": 0.0, "max": 0.0, "mean": 0.0}
    return {
        "min": float(min(lengths)),
        "max": float(max(lengths)),
        "mean": float(sum(lengths) / len(lengths)),
    }
