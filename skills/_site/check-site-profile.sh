#!/usr/bin/env bash
# Source the site profile named by HVANTK_SITE_PROFILE and report every required
# value that is missing. Exit 1 if any is. Used by every template before it submits.
#
#   check-site-profile.sh [--allow-missing-image]
#
# --allow-missing-image still requires HVANTK_SIF to be set, but does not require the
# image file to exist yet: use it before the first image build or promote, when no
# image has been produced on this site.
set -euo pipefail

ALLOW_MISSING_IMAGE=0
while (( $# )); do
  case "$1" in
    --allow-missing-image) ALLOW_MISSING_IMAGE=1; shift ;;
    *) echo "ERROR: unknown option $1" >&2; exit 2 ;;
  esac
done

: "${HVANTK_SITE_PROFILE:?export HVANTK_SITE_PROFILE=/path/to/site-profile.local.sh (copy skills/_site/site-profile.example.sh)}"
[[ -r "$HVANTK_SITE_PROFILE" ]] || { echo "ERROR: site profile not readable: $HVANTK_SITE_PROFILE" >&2; exit 1; }
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"

required=(HVANTK_REPO HVANTK_SIF HVANTK_CPU_PARTITION HVANTK_NODE_CPUS HVANTK_NODE_MEM_GB HVANTK_WORK)
missing=()
for v in "${required[@]}"; do
  [[ -n "${!v:-}" ]] || missing+=("$v")
done
if (( ${#missing[@]} )); then
  echo "ERROR: site profile $HVANTK_SITE_PROFILE leaves required values empty:" >&2
  printf '  %s\n' "${missing[@]}" >&2
  echo "Fill them in (see skills/_site/SITE-PROFILE.md); do not guess them." >&2
  exit 1
fi
[[ -d "$HVANTK_REPO" ]] || { echo "ERROR: HVANTK_REPO is not a directory: $HVANTK_REPO" >&2; exit 1; }
if [[ -e "$HVANTK_SIF" ]]; then
  :
elif (( ALLOW_MISSING_IMAGE )); then
  echo "NOTE: HVANTK_SIF does not exist yet: $HVANTK_SIF (run hvantk:hpc-container first)" >&2
else
  echo "ERROR: HVANTK_SIF not found: $HVANTK_SIF" >&2; exit 1
fi
[[ -d "$HVANTK_WORK" ]] || { echo "ERROR: HVANTK_WORK is not a directory: $HVANTK_WORK" >&2; exit 1; }
[[ "$HVANTK_NODE_CPUS" =~ ^[0-9]+$ && "$HVANTK_NODE_MEM_GB" =~ ^[0-9]+$ ]] \
  || { echo "ERROR: HVANTK_NODE_CPUS and HVANTK_NODE_MEM_GB must be integers" >&2; exit 1; }
echo "site profile OK: $HVANTK_SITE_PROFILE"
