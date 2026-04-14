"""Tests for the local corpus smoke-release builder."""

import json
from pathlib import Path
from types import SimpleNamespace

import ferritin
import ferritin.corpus_smoke as corpus_smoke


def _fake_structure(name: str):
    residues = [
        SimpleNamespace(name="GLY", serial_number=1, is_amino_acid=True, atoms=[]),
        SimpleNamespace(name="SER", serial_number=2, is_amino_acid=True, atoms=[]),
    ]
    chain = SimpleNamespace(id="A", residues=residues)
    return SimpleNamespace(identifier=name, chain_count=1, chains=[chain], residue_count=2, atom_count=10)


def test_build_local_corpus_smoke_release_orchestrates_pipeline(tmp_path, monkeypatch):
    loaded = [(0, _fake_structure("one")), (1, _fake_structure("two"))]

    def fake_batch_load_tolerant(paths, n_threads=None):
        return loaded

    def fake_batch_prepare(structures, n_threads=None):
        return [ferritin.PrepReport(hydrogens_added=1, converged=True) for _ in structures]

    def fake_build_structure_supervision_dataset_from_prepared(structures, prep_reports, out_dir, **kwargs):
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "prepared_structures.jsonl").write_text('{"record_id":"one"}\n{"record_id":"two"}\n', encoding="utf-8")
        sup = out / "supervision_release"
        sup.mkdir(parents=True, exist_ok=True)
        (sup / "release_manifest.json").write_text(json.dumps({"count_examples": 2, "lengths": {"mean": 2.0}}), encoding="utf-8")
        (sup / "failures.jsonl").write_text("", encoding="utf-8")
        (sup / "examples").mkdir(exist_ok=True)
        import numpy as np
        np.savez_compressed(
            sup / "examples" / "tensors.npz",
            seq_mask=np.ones((2, 2), dtype=np.float32),
            rigidgroups_gt_exists=np.ones((2, 2, 8), dtype=np.float32),
            pseudo_beta_mask=np.ones((2, 2), dtype=np.float32),
            chi_mask=np.ones((2, 2, 4), dtype=np.float32),
        )
        return out

    def fake_build_sequence_dataset(structures, out_dir, **kwargs):
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "release_manifest.json").write_text(json.dumps({"count_examples": 2, "lengths": {"mean": 2.0}}), encoding="utf-8")
        (out / "failures.jsonl").write_text("", encoding="utf-8")
        (out / "examples").mkdir(exist_ok=True)
        return out

    def fake_build_training_release(sequence_release_dir, structure_release_dir, out_dir, **kwargs):
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "release_manifest.json").write_text(json.dumps({"count_examples": 2, "split_counts": {"train": 1, "val": 1}}), encoding="utf-8")
        (out / "training_examples.jsonl").write_text(
            '{"record_id":"one","split":"train"}\n{"record_id":"two","split":"val"}\n',
            encoding="utf-8",
        )
        return out

    monkeypatch.setattr(corpus_smoke, "batch_load_tolerant", fake_batch_load_tolerant)
    monkeypatch.setattr(corpus_smoke, "batch_prepare", fake_batch_prepare)
    monkeypatch.setattr(corpus_smoke, "build_structure_supervision_dataset_from_prepared", fake_build_structure_supervision_dataset_from_prepared)
    monkeypatch.setattr(corpus_smoke, "build_sequence_dataset", fake_build_sequence_dataset)
    monkeypatch.setattr(corpus_smoke, "build_training_release", fake_build_training_release)

    root = ferritin.build_local_corpus_smoke_release(
        [tmp_path / "one.pdb", tmp_path / "two.pdb"],
        tmp_path / "smoke",
        release_id="smoke-v0",
        overwrite=True,
    )

    assert (root / "prepared" / "prepared_structures.jsonl").exists()
    assert (root / "sequence" / "release_manifest.json").exists()
    assert (root / "training" / "release_manifest.json").exists()
    assert (root / "corpus" / "corpus_release_manifest.json").exists()
    assert (root / "corpus" / "validation_report.json").exists()
