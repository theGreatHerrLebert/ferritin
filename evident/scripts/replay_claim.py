"""Replay a single claim by id from inside the EVIDENT image.

Looks up the claim's ``evidence.command`` from the manifest and runs
it as a subprocess in the source directory the claim points to. The
command resolves relative to ``source`` (default ``..``) under the
claim file's directory, so commands like ``pytest tests/oracle/...``
or ``python validation/charmm19_eef1_ball_oracle.py`` find their
dependencies under ``/workspace/proteon``.

Usage::

    python evident/scripts/replay_claim.py <claim-id> [-- extra args]

The ``-- extra args`` form forwards any extra tokens to the underlying
command, joined with the claim's ``evidence.command``. Useful for
``-k <pattern>`` or ``--maxfail=1`` style overrides at replay time.

Exit codes:
  0  command succeeded
  1  command failed (non-zero exit from the replayed command)
  64 usage / argument error
  65 claim id not found
  66 claim has no ``evidence.command``
"""
from __future__ import annotations

import argparse
import os
import pathlib
import shlex
import subprocess
import sys

import yaml

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
MANIFEST = REPO_ROOT / "evident" / "evident.yaml"


def _load_yaml(path: pathlib.Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def _all_claim_files() -> list[pathlib.Path]:
    """Return every YAML file under evident/claims/."""
    claims_dir = REPO_ROOT / "evident" / "claims"
    if not claims_dir.is_dir():
        return []
    return sorted(claims_dir.glob("*.yaml"))


def _find_claim(claim_id: str) -> tuple[pathlib.Path, dict] | None:
    for path in _all_claim_files():
        try:
            doc = _load_yaml(path)
        except yaml.YAMLError as e:
            print(
                f"warning: skipping unparseable claim file {path.name}: "
                f"{str(e).splitlines()[0]}",
                file=sys.stderr,
            )
            continue
        for claim in doc.get("claims") or []:
            if claim.get("id") == claim_id:
                return path, claim
    return None


def _resolve_source_dir(claim_yaml: pathlib.Path, claim: dict) -> pathlib.Path:
    """The directory the claim's evidence.command should run in.

    Defaults to ``..`` relative to the manifest (i.e. the repo root)
    when ``source`` is missing or set to ``..``.
    """
    source = claim.get("source", "..")
    base = (REPO_ROOT / "evident").resolve()
    return (base / source).resolve()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Replay a single EVIDENT claim by id."
    )
    parser.add_argument("claim_id", help="Claim id (matches `id` in the YAML).")
    parser.add_argument(
        "extra",
        nargs=argparse.REMAINDER,
        help="Optional extra tokens forwarded to the replayed command (after `--`).",
    )
    args = parser.parse_args()

    if not MANIFEST.exists():
        print(f"manifest not found: {MANIFEST}", file=sys.stderr)
        return 64

    found = _find_claim(args.claim_id)
    if found is None:
        print(f"claim not found: {args.claim_id}", file=sys.stderr)
        return 65
    claim_yaml, claim = found

    cmd = (claim.get("evidence") or {}).get("command")
    if not cmd:
        print(
            f"claim {args.claim_id}: no evidence.command — nothing to replay",
            file=sys.stderr,
        )
        return 66

    cwd = _resolve_source_dir(claim_yaml, claim)
    if not cwd.is_dir():
        print(f"source dir does not exist: {cwd}", file=sys.stderr)
        return 64

    extra = list(args.extra or [])
    if extra and extra[0] == "--":
        extra = extra[1:]

    cmdline = shlex.split(str(cmd)) + extra

    print(f"# replay claim: {args.claim_id}")
    print(f"# source:       {claim_yaml.relative_to(REPO_ROOT)}")
    print(f"# cwd:          {cwd.relative_to(REPO_ROOT)}")
    print(f"# command:      {' '.join(shlex.quote(a) for a in cmdline)}")
    sys.stdout.flush()

    try:
        result = subprocess.run(cmdline, cwd=cwd, env=os.environ.copy())
    except FileNotFoundError as e:
        print(f"command not found: {e}", file=sys.stderr)
        return 64
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
