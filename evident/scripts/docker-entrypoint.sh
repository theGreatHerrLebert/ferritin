#!/usr/bin/env bash
# EVIDENT image dispatch.
#
# Subcommands:
#   replay <claim-id> [-- extra args]   Run the evidence.command for the
#                                        named claim, in /workspace/proteon.
#   render <release-tag> [-- extra args] Freeze the claim bundle at tag time.
#                                        Wraps evident/scripts/lock_release_replays.py.
#   list | validate | <other>            Pass-through to the vendored
#                                        `evident` CLI for backwards compat.
#
# Examples:
#   docker run --rm <image> list evident/evident.yaml
#   docker run --rm <image> replay proteon-charmm19-vs-ball-corpus-1k-pdbs
#   docker run --rm -v "$(pwd)/out:/workspace/proteon/evident/reports/v1" \
#               <image> render v1
#
# The image is intentionally a thin envelope: every subcommand resolves
# back to scripts/files that already live in the repo, so a third party
# can replicate the same run without docker simply by checking out the
# tag and invoking the same Python entrypoint.

set -euo pipefail

cd /workspace/proteon

# v0.2.0 data-mount contract: image bundles tools + executors, user mounts
# data via --bind. Auto-export the well-known env vars when the convention
# mounts exist so that `docker run -v ./pdbs:/data/pdbs -v ./out:/data/out
# <image> replay <claim-id>` works with no -e flags. User-supplied
# PROTEON_CORPUS_DIR / PROTEON_OUTPUT_DIR override the defaults below.
: "${PROTEON_CORPUS_DIR:=/data/pdbs}"
: "${PROTEON_OUTPUT_DIR:=/data/out}"
if [[ -d "$PROTEON_CORPUS_DIR" ]]; then export PROTEON_CORPUS_DIR; else unset PROTEON_CORPUS_DIR; fi
if [[ -d "$PROTEON_OUTPUT_DIR" ]]; then export PROTEON_OUTPUT_DIR; else unset PROTEON_OUTPUT_DIR; fi
if [[ -x /usr/local/bin/USalign && -z "${USALIGN_BIN:-}" ]]; then export USALIGN_BIN=/usr/local/bin/USalign; fi

usage() {
    cat <<USAGE
Usage:
  $0 replay <claim-id> [-- extra args]
  $0 render <release-tag> [-- extra args]
  $0 list|validate|... (forwarded to the evident CLI)
USAGE
}

if [[ $# -lt 1 ]]; then
    exec /usr/local/bin/evident list evident/evident.yaml
fi

case "$1" in
    replay)
        shift
        if [[ $# -lt 1 ]]; then
            echo "replay: missing <claim-id>" >&2
            usage >&2
            exit 64
        fi
        exec python3 evident/scripts/replay_claim.py "$@"
        ;;
    render)
        shift
        if [[ $# -lt 1 ]]; then
            echo "render: missing <release-tag>" >&2
            usage >&2
            exit 64
        fi
        tag="$1"
        shift
        exec python3 evident/scripts/lock_release_replays.py --release "$tag" "$@"
        ;;
    -h|--help|help)
        usage
        exit 0
        ;;
    *)
        exec /usr/local/bin/evident "$@"
        ;;
esac
