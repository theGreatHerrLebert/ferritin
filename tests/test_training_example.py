"""Tests for joined training example artifacts."""

import json
from types import SimpleNamespace

import ferritin


def _atom(name, xyz):
    return SimpleNamespace(name=name, pos=tuple(float(x) for x in xyz))


def _fake_structure(chain_id="A"):
    residues = [
        SimpleNamespace(
            name="GLY",
            serial_number=1,
            is_amino_acid=True,
            atoms=[_atom("N", (0, 0, 0)), _atom("CA", (1, 0, 0)), _atom("C", (1.8, 1, 0)), _atom("O", (1.8, 2.1, 0))],
        ),
        SimpleNamespace(
            name="SER",
            serial_number=2,
            is_amino_acid=True,
            atoms=[
                _atom("N", (2.6, 0.8, 0.8)),
                _atom("CA", (3.5, 1.6, 1.0)),
                _atom("C", (4.7, 0.9, 1.4)),
                _atom("O", (5.0, -0.2, 1.1)),
                _atom("CB", (3.2, 2.9, 1.8)),
                _atom("OG", (2.1, 3.5, 1.4)),
            ],
        ),
    ]
    chain = SimpleNamespace(id=chain_id, residues=residues)
    return SimpleNamespace(identifier="fake", chain_count=1, chains=[chain])


class TestTrainingExample:
    def test_join_training_release_from_sequence_and_structure(self, tmp_path):
        structures = [_fake_structure("A")]
        seq_release = ferritin.build_sequence_dataset(
            structures,
            tmp_path / "seq_release",
            release_id="seq-v0",
            record_ids=["fake:A"],
        )
        prep = ferritin.PrepReport(hydrogens_added=2, converged=True)
        struc_release_root = ferritin.build_structure_supervision_dataset_from_prepared(
            structures,
            [prep],
            tmp_path / "struc_release_root",
            release_id="struc-v0",
            record_ids=["fake:A"],
        )
        train_release = ferritin.build_training_release(
            seq_release,
            struc_release_root / "supervision_release",
            tmp_path / "train_release",
            release_id="train-v0",
            split_assignments={"fake:A": "train"},
            crop_metadata={"fake:A": (0, 2)},
            weights={"fake:A": 0.5},
        )

        manifest = json.loads((train_release / "release_manifest.json").read_text(encoding="utf-8"))
        assert manifest["release_id"] == "train-v0"
        assert manifest["count_examples"] == 1
        assert manifest["split_counts"]["train"] == 1

        row = json.loads((train_release / "training_examples.jsonl").read_text(encoding="utf-8").strip())
        assert row["record_id"] == "fake:A"
        assert row["split"] == "train"
        assert row["crop_start"] == 0
        assert row["crop_stop"] == 2
        assert row["weight"] == 0.5
