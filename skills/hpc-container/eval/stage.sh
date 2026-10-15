#!/usr/bin/env bash
# Stage for the hpc-container eval: nothing to download. Confirms an unprivileged build
# method exists on this site and creates the cache directories under the work root.
#
#   HVANTK_SITE_PROFILE=... bash skills/hpc-container/eval/stage.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/hpc-container"
mkdir -p "$STAGE" "$HVANTK_WORK/containers/cache" "$HVANTK_WORK/containers/tmp"
if [[ "$(grep -c "^$USER:" /etc/subuid 2>/dev/null || true)" -gt 0 ]]; then
  echo "build method: fakeroot (subuid entry present)"
elif command -v proot >/dev/null 2>&1; then
  echo "build method: unprivileged + proot ($(command -v proot))"
else
  echo "ERROR: no subuid entry and no proot on PATH; see skills/hpc-container/references/build-troubleshooting.md" >&2
  exit 1
fi
: > "$STAGE/staged.sha256"   # no staged input; the file exists for the contract
echo "staged under $STAGE"
