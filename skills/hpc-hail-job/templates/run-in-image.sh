#!/usr/bin/env bash
# Run one command inside the hvantk container image on the current node, with the
# checkout on PYTHONPATH (so the checkout's code runs, and the files it writes, such as
# a drift fingerprint or the rebuild ledger, land in the checkout), and the work root,
# the checkout and the eval data root bound. For login-node work: staging downloads,
# drift probes, bookkeeping, small verifications. Anything Hail-heavy goes through
# submit-hail-job.sh as a job.
#
#   run-in-image.sh [--sif <image>] -- <hvantk ...|python ...>
#
# Reads every site value from $HVANTK_SITE_PROFILE (skills/_site/SITE-PROFILE.md).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SITE_DIR="$(cd "$HERE/../../_site" && pwd -P)"
bash "$SITE_DIR/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"

SIF="$HVANTK_SIF"
while (( $# )); do
  case "$1" in
    --sif) SIF="$2"; shift 2 ;;
    --) shift; break ;;
    *) echo "ERROR: unknown option $1 (put the command after --)" >&2; exit 2 ;;
  esac
done
(( $# > 0 )) || { echo "ERROR: no command after --" >&2; exit 2; }
case "$1" in hvantk|python|python3) ;; *) echo "ERROR: command must start with hvantk or python, got '$1'" >&2; exit 2 ;; esac
[[ -e "$SIF" ]] || { echo "ERROR: image not found: $SIF" >&2; exit 2; }
CONTAINER_BIN="${HVANTK_CONTAINER_BIN:-singularity}"
command -v "$CONTAINER_BIN" >/dev/null || { echo "ERROR: $CONTAINER_BIN not on PATH" >&2; exit 2; }
[[ -d "$HVANTK_REPO/hvantk" ]] || { echo "ERROR: checkout $HVANTK_REPO has no hvantk/ package; the image's own copy would run instead" >&2; exit 2; }

BINDS=(-B "$HVANTK_WORK:$HVANTK_WORK")
bind_outside_work() {
  local p="$1"
  [[ -d "$p" ]] || return 0
  case "$p" in "$HVANTK_WORK"|"$HVANTK_WORK"/*) ;; *) BINDS+=(-B "$p:$p") ;; esac
}
bind_outside_work "$HVANTK_REPO"
[[ -n "${HVANTK_EVAL_DATA:-}" ]] && bind_outside_work "$HVANTK_EVAL_DATA"
[[ -n "${HVANTK_BIND:-}" ]] && BINDS+=(-B "$HVANTK_BIND")

# The checkout's commit, for the builder_commit a build records: the image has no git.
# Same rule as hvantk's own resolver: HEAD, only when HVANTK_REPO is the top of its own
# git checkout, with -dirty for uncommitted edits under hvantk/. Empty when git cannot
# answer here; the build then records the hvantk version.
COMMIT=""
if command -v git >/dev/null 2>&1 \
  && [[ "$(git -C "$HVANTK_REPO" rev-parse --show-toplevel 2>/dev/null)" -ef "$HVANTK_REPO" ]]; then
  COMMIT="$(git -C "$HVANTK_REPO" rev-parse HEAD 2>/dev/null)" || COMMIT=""
  if [[ -n "$COMMIT" && -n "$(git -C "$HVANTK_REPO" status --porcelain --untracked-files=no -- hvantk 2>/dev/null)" ]]; then
    COMMIT="$COMMIT-dirty"
  fi
fi

exec "$CONTAINER_BIN" exec "${BINDS[@]}" \
  --env PYTHONPATH="$HVANTK_REPO" --env HVANTK_BUILDER_COMMIT="$COMMIT" \
  --env JAVA_TOOL_OPTIONS="-Duser.language=en -Duser.country=US" \
  "$SIF" "$@"
