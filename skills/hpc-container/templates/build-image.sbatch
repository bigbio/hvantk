#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=hvantk-build-image
# Resources, partition, time, output and mail are passed on the command line by
# submit-build-image.sh, which also exports the git facts (compute nodes often have no
# git). No site value is written here.
#
# Builds the hvantk image from the checkout's def file and lock file, without root
# (fakeroot where /etc/subuid allows it, a static proot otherwise), gates the Java and
# Hail versions inside the image, runs `hvantk utils check-install` through node-local
# scratch, and writes <image>.provenance.json next to the image before the final gate.
set -euo pipefail

: "${HVANTK_SITE_PROFILE:?HVANTK_SITE_PROFILE not set}"
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
: "${HVANTK_REPO:?}" "${HVANTK_WORK:?}"
: "${BUILD_OUT:?BUILD_OUT not set (submit-build-image.sh exports it)}"
: "${BUILD_DEF:?}" "${BUILD_GIT_SHA:?}" "${BUILD_GIT_BRANCH:?}" "${BUILD_GIT_DIRTY:?}"

CONTAINER_BIN="${HVANTK_CONTAINER_BIN:-singularity}"
command -v "$CONTAINER_BIN" >/dev/null || { echo "ERROR: $CONTAINER_BIN not on PATH" >&2; exit 2; }
[[ -f "$BUILD_DEF" ]] || { echo "ERROR: def file not found: $BUILD_DEF" >&2; exit 2; }
[[ -f "$HVANTK_REPO/poetry.lock" ]] || { echo "ERROR: no poetry.lock in $HVANTK_REPO" >&2; exit 2; }
mkdir -p "$(dirname "$BUILD_OUT")"

# --- cache and temporary files under the work root, never the home quota ---
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-$HVANTK_WORK/containers/cache}"
export SINGULARITY_TMPDIR="${SINGULARITY_TMPDIR:-$HVANTK_WORK/containers/tmp}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$SINGULARITY_CACHEDIR}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$SINGULARITY_TMPDIR}"
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

# --- unprivileged build method ---
BUILD_FLAGS=()
if [[ "$(grep -c "^$USER:" /etc/subuid 2>/dev/null || true)" -gt 0 ]]; then
  METHOD="fakeroot"; BUILD_FLAGS=(--fakeroot)
elif command -v proot >/dev/null 2>&1; then
  METHOD="unprivileged + proot"
else
  echo "ERROR: neither a /etc/subuid entry (fakeroot) nor a proot on PATH; see references/build-troubleshooting.md" >&2
  exit 2
fi

# --- egress: the base image and the wheels are pulled during the build ---
CODE="$(curl -sS --max-time 20 -o /dev/null -w '%{http_code}' https://registry-1.docker.io/v2/ || true)"
[[ -n "$CODE" && "$CODE" != "000" ]] || { echo "ERROR: no route to the container registry from $(hostname); build where egress exists" >&2; exit 2; }

LOCK_MD5="$(md5sum "$HVANTK_REPO/poetry.lock" | cut -d' ' -f1)"
echo "== hvantk image build =="
echo "date:    $(date -Is)"
echo "host:    $(hostname)  job=${SLURM_JOB_ID:-none}"
echo "repo:    $HVANTK_REPO @ $BUILD_GIT_SHA ($BUILD_GIT_BRANCH, $BUILD_GIT_DIRTY)"
echo "def:     $BUILD_DEF"
echo "out:     $BUILD_OUT"
echo "method:  $METHOD  tool: $("$CONTAINER_BIN" --version)"
echo "lock:    md5 $LOCK_MD5"
echo "========================"

# %files paths in the def are relative to the build directory: the repository root.
cd "$HVANTK_REPO"
"$CONTAINER_BIN" build -F "${BUILD_FLAGS[@]}" "$BUILD_OUT" "$BUILD_DEF"
[[ -f "$BUILD_OUT" ]] || { echo "ERROR: build reported success but $BUILD_OUT is missing" >&2; exit 1; }
du -h "$BUILD_OUT"

# --- gates inside the image ---
JAVA_VERSION="$("$CONTAINER_BIN" exec "$BUILD_OUT" java -version 2>&1 | head -1 | sed -E 's/.*"([^"]+)".*/\1/')"
[[ "$JAVA_VERSION" =~ ^(1\.8|11)\. ]] || { echo "ERROR: Java inside the image is $JAVA_VERSION; Hail needs 8 or 11" >&2; exit 1; }
HAIL_VERSION="$("$CONTAINER_BIN" exec "$BUILD_OUT" python3 -c 'import hail; print(hail.__version__)')"
echo "java:    $JAVA_VERSION"
echo "hail:    $HAIL_VERSION"

# --- install check through node-local scratch (Hail cannot start without it) ---
SCRATCH="${HVANTK_SCRATCH_ROOT:-${SLURM_TMPDIR:-/tmp}}/hvantk-build-$USER-${SLURM_JOB_ID:-$$}"
mkdir -p "$SCRATCH"
set +e
"$CONTAINER_BIN" exec -B "$SCRATCH:$SCRATCH" \
  --env SPARK_LOCAL_DIRS="$SCRATCH" --env TMPDIR="$SCRATCH" --env HAIL_TMPDIR="$SCRATCH" \
  "$BUILD_OUT" hvantk utils check-install
CHECK_RC=$?
set -e
rm -rf "$SCRATCH"

# --- provenance record, written before the final gate so a failure is recorded too ---
cat > "${BUILD_OUT}.provenance.json" <<JSON
{
  "sif": "$BUILD_OUT",
  "built_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "slurm_job": "${SLURM_JOB_ID:-}",
  "host": "$(hostname)",
  "repo": "$HVANTK_REPO",
  "git_commit": "$BUILD_GIT_SHA",
  "git_branch": "$BUILD_GIT_BRANCH",
  "worktree": "$BUILD_GIT_DIRTY",
  "poetry_lock_md5": "$LOCK_MD5",
  "container_tool": "$("$CONTAINER_BIN" --version)",
  "build_method": "$METHOD",
  "hail_version": "$HAIL_VERSION",
  "java_version": "$JAVA_VERSION",
  "check_install_rc": $CHECK_RC
}
JSON
echo "record:  ${BUILD_OUT}.provenance.json"
[[ "$CHECK_RC" -eq 0 ]] || { echo "ERROR: hvantk utils check-install failed inside the new image (rc=$CHECK_RC)" >&2; exit 1; }
echo "DONE: $BUILD_OUT"
