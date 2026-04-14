"""Smoke pipeline for building a small real-file corpus release."""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

from .corpus_release import build_corpus_release_manifest
from .corpus_validation import validate_corpus_release
from .failure_taxonomy import PARSE_ERROR
from .io import batch_load_tolerant
from .prepare import batch_prepare
from .sequence_release import build_sequence_dataset
from .supervision_dataset import build_structure_supervision_dataset_from_prepared
from .supervision_release import FailureRecord
from .training_example import build_training_release


def build_local_corpus_smoke_release(
    paths: Sequence[str | Path],
    out_dir: str | Path,
    *,
    release_id: str,
    code_rev: Optional[str] = None,
    config_rev: Optional[str] = None,
    prep_policy_version: Optional[str] = None,
    split_policy_version: Optional[str] = None,
    n_threads: Optional[int] = None,
    overwrite: bool = False,
) -> Path:
    """Build a small end-to-end corpus release from local structure files.

    This is the intended smoke path for validating the full data-release stack
    on real PDB/mmCIF inputs already available on disk.
    """
    root = Path(out_dir)
    if root.exists() and not overwrite:
        raise FileExistsError(f"{root} already exists")
    root.mkdir(parents=True, exist_ok=True)

    path_list = [Path(p) for p in paths]
    loaded_pairs = batch_load_tolerant(path_list, n_threads=n_threads)
    loaded_indices = [idx for idx, _ in loaded_pairs]
    loaded_structures = [structure for _, structure in loaded_pairs]
    loaded_paths = [path_list[idx] for idx in loaded_indices]
    record_ids = [path.stem for path in loaded_paths]
    source_ids = [str(path) for path in loaded_paths]

    # Capture ingestion failures. `batch_load_tolerant` silently drops
    # bad paths; if we only forward what it returned we lose visibility
    # into partial ingestion. Roadmap Section 7 requires every failure
    # mode to be machine-readable, so we write one FailureRecord per
    # dropped path and merge the breakdown into the top-level manifest.
    loaded_set = set(loaded_indices)
    dropped_indices = [i for i in range(len(path_list)) if i not in loaded_set]
    ingestion_failures_path = _write_ingestion_failures(
        root,
        [path_list[i] for i in dropped_indices],
        code_rev=code_rev,
        config_rev=config_rev,
    )

    prep_reports = batch_prepare(loaded_structures, n_threads=n_threads)

    prepared_root = build_structure_supervision_dataset_from_prepared(
        loaded_structures,
        prep_reports,
        root / "prepared",
        release_id=f"{release_id}-structure",
        record_ids=record_ids,
        source_ids=source_ids,
        code_rev=code_rev,
        config_rev=config_rev,
        provenance={"input_paths": [str(p) for p in loaded_paths]},
        overwrite=True,
    )
    sequence_root = build_sequence_dataset(
        loaded_structures,
        root / "sequence",
        release_id=f"{release_id}-sequence",
        record_ids=record_ids,
        source_ids=source_ids,
        code_rev=code_rev,
        config_rev=config_rev,
        provenance={"input_paths": [str(p) for p in loaded_paths]},
        overwrite=True,
    )
    split_assignments = _default_split_assignments(record_ids)
    training_root = build_training_release(
        sequence_root,
        prepared_root / "supervision_release",
        root / "training",
        release_id=f"{release_id}-training",
        split_assignments=split_assignments,
        code_rev=code_rev,
        config_rev=config_rev,
        provenance={"input_paths": [str(p) for p in loaded_paths]},
        overwrite=True,
    )
    corpus_root = build_corpus_release_manifest(
        root / "corpus",
        release_id=release_id,
        prepared_manifest=prepared_root / "prepared_structures.jsonl",
        sequence_release=sequence_root,
        structure_release=prepared_root / "supervision_release",
        training_release=training_root,
        ingestion_failures=ingestion_failures_path,
        code_rev=code_rev,
        config_rev=config_rev,
        prep_policy_version=prep_policy_version,
        split_policy_version=split_policy_version,
        provenance={
            "input_paths": [str(p) for p in path_list],
            "loaded_paths": [str(p) for p in loaded_paths],
            "dropped_paths": [str(path_list[i]) for i in dropped_indices],
        },
        overwrite=True,
    )
    validate_corpus_release(
        corpus_root / "corpus_release_manifest.json",
        out_path=corpus_root / "validation_report.json",
    )
    return root


def _write_ingestion_failures(
    root: Path,
    dropped_paths: List[Path],
    *,
    code_rev: Optional[str],
    config_rev: Optional[str],
) -> Optional[Path]:
    """Write `FailureRecord` JSONL for inputs that didn't load.

    Returns `None` when nothing was dropped — callers / the corpus
    manifest builder treat `None` as "no ingestion failures to merge",
    so we don't produce empty files.
    """
    if not dropped_paths:
        return None
    ingestion_path = root / "ingestion_failures.jsonl"
    with ingestion_path.open("w", encoding="utf-8") as handle:
        for path in dropped_paths:
            record = FailureRecord(
                record_id=path.stem,
                stage="raw_intake",
                failure_class=PARSE_ERROR,
                message=(
                    f"batch_load_tolerant could not parse {path}; "
                    "the upstream loader silently dropped this input"
                ),
                source_id=str(path),
                code_rev=code_rev,
                config_rev=config_rev,
                provenance={"stage_entry_point": "corpus_smoke.batch_load_tolerant"},
            )
            handle.write(json.dumps(asdict(record), separators=(",", ":")))
            handle.write("\n")
    return ingestion_path


def _default_split_assignments(record_ids: Iterable[str]) -> dict[str, str]:
    ordered = list(record_ids)
    if not ordered:
        return {}
    if len(ordered) == 1:
        return {ordered[0]: "train"}
    assignments = {record_id: "train" for record_id in ordered}
    assignments[ordered[-1]] = "val"
    return assignments
