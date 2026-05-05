"""Guard the v0.2.0 EVIDENT data-mount contract on every release-tier runner.

Each release-tier oracle runner under ``validation/`` must honour
``PROTEON_CORPUS_DIR`` (input directory) and ``PROTEON_OUTPUT_DIR``
(output directory) env vars when set, falling back to its historical
default when unset. This file imports each runner module with controlled
env vars and asserts the resulting module-level path constants.

This is a pure path-resolution check — no oracle calls, no heavy compute.
It runs in milliseconds and exists to catch the v0.1.x failure mode where
runners hardcoded monster3-only paths and broke under the v0.2.0 image's
``--bind`` model.

Heavy oracle imports (openmm, pdbfixer, ball-py, pebble) trigger when a
runner module loads. Tests skip gracefully if any required dep is missing,
so the contract check still runs in lighter test envs.
"""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATION = REPO_ROOT / "validation"


def _load_runner(module_path: Path):
    """Import a runner .py from disk as a fresh module.

    Bypasses the import cache so each call re-evaluates the module-level
    env-var lookups against the current ``os.environ``.
    """
    spec = importlib.util.spec_from_file_location(
        f"_contract_test_{module_path.stem}", module_path
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _clear_proteon_env(monkeypatch):
    for key in list(os.environ):
        if key.startswith("PROTEON_") or key == "USALIGN_BIN":
            monkeypatch.delenv(key, raising=False)


# Per-runner contract: (path, expected output filename within PROTEON_OUTPUT_DIR,
# required-to-import dep, dep-name-for-skip-message)
FOLD_PRES_RUNNERS = [
    ("tm_fold_preservation.py",          "tm_fold_preservation.jsonl",          "proteon"),
    ("tm_fold_preservation_amber.py",    "tm_fold_preservation_amber.jsonl",    "proteon"),
    ("tm_fold_preservation_openmm.py",   "tm_fold_preservation_openmm.jsonl",   "openmm"),
    ("tm_fold_preservation_openmm_amber.py", "tm_fold_preservation_openmm_amber.jsonl", "openmm"),
]


@pytest.mark.parametrize("runner_file,canonical_out,required_dep", FOLD_PRES_RUNNERS)
def test_fold_preservation_runner_honors_env(
    monkeypatch, tmp_path, runner_file, canonical_out, required_dep
):
    pytest.importorskip(required_dep)
    _clear_proteon_env(monkeypatch)
    pdbs = tmp_path / "pdbs"
    pdbs.mkdir()
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setenv("PROTEON_CORPUS_DIR", str(pdbs))
    monkeypatch.setenv("PROTEON_OUTPUT_DIR", str(out))

    m = _load_runner(VALIDATION / runner_file)

    assert Path(m.PDB_DIR) == pdbs, f"{runner_file}: PDB_DIR ignores PROTEON_CORPUS_DIR"
    assert Path(m.OUT) == out / canonical_out, (
        f"{runner_file}: OUT does not resolve to PROTEON_OUTPUT_DIR/{canonical_out}"
    )


@pytest.mark.parametrize("runner_file,canonical_out,required_dep", FOLD_PRES_RUNNERS)
def test_fold_preservation_runner_default_unchanged(
    monkeypatch, runner_file, canonical_out, required_dep
):
    """Without env vars set, runners must keep their historical defaults so
    monster3 batch invocations and the legacy scp-then-join workflow keep
    working unchanged."""
    pytest.importorskip(required_dep)
    _clear_proteon_env(monkeypatch)

    m = _load_runner(VALIDATION / runner_file)

    expected_pdb = Path("/globalscratch/dateschn/proteon-benchmark/pdbs_50k")
    expected_out = Path("/globalscratch/dateschn/proteon-benchmark") / canonical_out
    assert Path(m.PDB_DIR) == expected_pdb, (
        f"{runner_file}: default PDB_DIR drifted from monster3 path"
    )
    assert Path(m.OUT) == expected_out, (
        f"{runner_file}: default OUT drifted from monster3 path"
    )


def test_join_fold_preservation_honors_env(monkeypatch, tmp_path):
    _clear_proteon_env(monkeypatch)
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setenv("PROTEON_OUTPUT_DIR", str(out))

    m = _load_runner(VALIDATION / "fold_preservation" / "join_fold_preservation.py")

    assert Path(m.SRC_DIR) == out, "join: SRC_DIR ignores PROTEON_OUTPUT_DIR"
    assert Path(m.DEST_DIR) == out, "join: DEST_DIR ignores PROTEON_OUTPUT_DIR"


def test_join_fold_preservation_default_unchanged(monkeypatch):
    _clear_proteon_env(monkeypatch)

    m = _load_runner(VALIDATION / "fold_preservation" / "join_fold_preservation.py")

    here = (VALIDATION / "fold_preservation").resolve()
    assert Path(m.SRC_DIR) == here, "join: default SRC_DIR drifted"
    assert Path(m.DEST_DIR) == here, "join: default DEST_DIR drifted"


def test_charmm_oracle_v02_synonyms(monkeypatch, tmp_path):
    """charmm19_eef1_ball_oracle.py accepts BOTH the legacy
    PROTEON_PDB_DIR / PROTEON_CHARMM_ORACLE_OUT and the new universal
    PROTEON_CORPUS_DIR / PROTEON_OUTPUT_DIR. Legacy first."""
    pytest.importorskip("proteon")
    pytest.importorskip("ball")
    pytest.importorskip("pebble")
    _clear_proteon_env(monkeypatch)
    pdbs = tmp_path / "pdbs"
    pdbs.mkdir()
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setenv("PROTEON_CORPUS_DIR", str(pdbs))
    monkeypatch.setenv("PROTEON_OUTPUT_DIR", str(out))

    m = _load_runner(VALIDATION / "charmm19_eef1_ball_oracle.py")

    assert Path(m.PDB_DIR) == pdbs
    assert Path(m.OUT) == out / "charmm19_eef1_ball_oracle.jsonl"


def test_charmm_oracle_legacy_env_wins(monkeypatch, tmp_path):
    """When BOTH legacy and v0.2.0 env vars are set, legacy takes precedence
    so existing monster3 batch scripts never silently change behaviour."""
    pytest.importorskip("proteon")
    pytest.importorskip("ball")
    pytest.importorskip("pebble")
    _clear_proteon_env(monkeypatch)
    legacy_dir = tmp_path / "legacy_pdbs"
    legacy_dir.mkdir()
    v02_dir = tmp_path / "v02_pdbs"
    v02_dir.mkdir()
    legacy_out = tmp_path / "legacy.jsonl"
    v02_out_dir = tmp_path / "v02_out"
    v02_out_dir.mkdir()

    monkeypatch.setenv("PROTEON_PDB_DIR", str(legacy_dir))
    monkeypatch.setenv("PROTEON_CORPUS_DIR", str(v02_dir))
    monkeypatch.setenv("PROTEON_CHARMM_ORACLE_OUT", str(legacy_out))
    monkeypatch.setenv("PROTEON_OUTPUT_DIR", str(v02_out_dir))

    m = _load_runner(VALIDATION / "charmm19_eef1_ball_oracle.py")

    assert Path(m.PDB_DIR) == legacy_dir
    assert Path(m.OUT) == legacy_out


def test_run_validation_usalign_env(monkeypatch, tmp_path):
    """run_validation.py honours USALIGN_BIN env var (image vendors USAlign
    at /usr/local/bin/USalign; source-tree dev keeps the
    /scratch/TMAlign/USAlign/ default)."""
    pytest.importorskip("proteon")
    _clear_proteon_env(monkeypatch)
    fake_bin = tmp_path / "USalign"
    fake_bin.write_text("#!/bin/sh\nexit 0\n")
    fake_bin.chmod(0o755)
    monkeypatch.setenv("USALIGN_BIN", str(fake_bin))

    m = _load_runner(VALIDATION / "run_validation.py")

    assert Path(m.USALIGN_BIN) == fake_bin
    assert m.HAS_USALIGN is True
