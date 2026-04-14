"""Framework-neutral joined training examples."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from .sequence_example import SequenceExample
from .sequence_export import load_sequence_examples
from .supervision import StructureSupervisionExample
from .supervision_export import load_structure_supervision_examples


@dataclass
class TrainingExample:
    """Thin join artifact over sequence and structure examples."""

    record_id: str
    source_id: Optional[str]
    chain_id: str
    split: str
    crop_start: Optional[int] = None
    crop_stop: Optional[int] = None
    weight: float = 1.0
    sequence: SequenceExample | None = None
    structure: StructureSupervisionExample | None = None


@dataclass
class TrainingReleaseManifest:
    """Shared manifest linking sequence and structure releases."""

    release_id: str
    artifact_type: str = "release_manifest"
    format: str = "ferritin.training_example.v0"
    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    code_rev: Optional[str] = None
    config_rev: Optional[str] = None
    sequence_release: str = ""
    structure_release: str = ""
    count_examples: int = 0
    split_counts: Dict[str, int] = field(default_factory=dict)
    examples_file: str = "training_examples.jsonl"
    provenance: Dict[str, object] = field(default_factory=dict)


def join_training_examples(
    sequence_examples: Sequence[SequenceExample],
    structure_examples: Sequence[StructureSupervisionExample],
    *,
    split_assignments: Optional[Dict[str, str]] = None,
    crop_metadata: Optional[Dict[str, tuple[int, int]]] = None,
    weights: Optional[Dict[str, float]] = None,
) -> List[TrainingExample]:
    """Join sequence and structure artifacts by `record_id`."""
    seq_by_id = {ex.record_id: ex for ex in sequence_examples}
    struc_by_id = {ex.record_id: ex for ex in structure_examples}
    shared_ids = sorted(set(seq_by_id).intersection(struc_by_id))

    out: List[TrainingExample] = []
    for record_id in shared_ids:
        seq = seq_by_id[record_id]
        struc = struc_by_id[record_id]
        split = (split_assignments or {}).get(record_id, "train")
        crop = (crop_metadata or {}).get(record_id)
        out.append(
            TrainingExample(
                record_id=record_id,
                source_id=seq.source_id or struc.source_id,
                chain_id=seq.chain_id,
                split=split,
                crop_start=None if crop is None else int(crop[0]),
                crop_stop=None if crop is None else int(crop[1]),
                weight=float((weights or {}).get(record_id, 1.0)),
                sequence=seq,
                structure=struc,
            )
        )
    return out


def build_training_release(
    sequence_release_dir: str | Path,
    structure_release_dir: str | Path,
    out_dir: str | Path,
    *,
    release_id: str,
    split_assignments: Optional[Dict[str, str]] = None,
    crop_metadata: Optional[Dict[str, tuple[int, int]]] = None,
    weights: Optional[Dict[str, float]] = None,
    code_rev: Optional[str] = None,
    config_rev: Optional[str] = None,
    provenance: Optional[Dict[str, object]] = None,
    overwrite: bool = False,
) -> Path:
    """Build a training release by joining sequence and structure releases."""
    sequence_examples = load_sequence_examples(Path(sequence_release_dir) / "examples")
    structure_examples = load_structure_supervision_examples(Path(structure_release_dir) / "examples")
    training_examples = join_training_examples(
        sequence_examples,
        structure_examples,
        split_assignments=split_assignments,
        crop_metadata=crop_metadata,
        weights=weights,
    )

    root = Path(out_dir)
    if root.exists() and not overwrite:
        raise FileExistsError(f"{root} already exists")
    root.mkdir(parents=True, exist_ok=True)

    rows = []
    split_counts: Dict[str, int] = {}
    for ex in training_examples:
        rows.append(
            {
                "record_id": ex.record_id,
                "source_id": ex.source_id,
                "chain_id": ex.chain_id,
                "split": ex.split,
                "crop_start": ex.crop_start,
                "crop_stop": ex.crop_stop,
                "weight": ex.weight,
            }
        )
        split_counts[ex.split] = split_counts.get(ex.split, 0) + 1

    with (root / "training_examples.jsonl").open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, separators=(",", ":")))
            handle.write("\n")

    manifest = TrainingReleaseManifest(
        release_id=release_id,
        code_rev=code_rev,
        config_rev=config_rev,
        sequence_release=str(Path(sequence_release_dir)),
        structure_release=str(Path(structure_release_dir)),
        count_examples=len(training_examples),
        split_counts=split_counts,
        provenance=dict(provenance or {}),
    )
    (root / "release_manifest.json").write_text(json.dumps(asdict(manifest), indent=2), encoding="utf-8")
    return root
