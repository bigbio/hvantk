#!/usr/bin/env bash
# Make a validated image the default one, keeping the old default under a dated name.
#
#   promote-image.sh <image.sif> [--default <path, default: the profile's HVANTK_SIF>] [--no-verify] [--name <job name>]
#
# Refuses while any of your jobs is running, and unless <image.sif>.provenance.json says
# check_install_rc is 0. Afterwards submits a verification job through the default name
# unless --no-verify is given. Login-node work.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" --allow-missing-image >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"

NEW="${1:?usage: promote-image.sh <image.sif> [--default <path>] [--no-verify]}"; shift
DEFAULT="$HVANTK_SIF"; VERIFY=1; NAME="promote-verify"
while (( $# )); do
  case "$1" in
    --default) DEFAULT="$2"; shift 2 ;;
    --no-verify) VERIFY=0; shift ;;
    --name) NAME="$2"; shift 2 ;;
    *) echo "ERROR: unknown option $1" >&2; exit 2 ;;
  esac
done
[[ -f "$NEW" ]] || { echo "ERROR: image not found: $NEW" >&2; exit 2; }
RECORD="$NEW.provenance.json"
[[ -f "$RECORD" ]] || { echo "ERROR: no provenance record next to $NEW; build it with build-image.sbatch" >&2; exit 2; }
grep -Eq '"check_install_rc": *0\b' "$RECORD" || { echo "ERROR: $RECORD does not record a passing install check; not promoting" >&2; exit 1; }
QUEUE="$(squeue -u "$USER" -h -o %i 2>/dev/null)" || { echo "ERROR: cannot query the queue; not promoting blind" >&2; exit 1; }
if [[ -n "$QUEUE" ]]; then
  echo "ERROR: you have jobs in the queue; a job may be using the default image. Wait for them." >&2; exit 1
fi
NEW_REAL="$(readlink -f "$NEW")"
[[ "$NEW_REAL" != "$(readlink -f "$DEFAULT" 2>/dev/null || true)" ]] || { echo "already the default: $NEW"; exit 0; }

if [[ -e "$DEFAULT" && ! -L "$DEFAULT" ]]; then
  # A real file as the default: keep it under a dated name with its record.
  KEEP="$(dirname "$DEFAULT")/hvantk-$(date -u +%Y-%m-%d).sif"
  [[ -e "$KEEP" ]] && KEEP="${KEEP%.sif}-$(date -u +%H%M%S).sif"
  mv "$DEFAULT" "$KEEP"
  [[ -f "$DEFAULT.provenance.json" ]] && mv "$DEFAULT.provenance.json" "$KEEP.provenance.json"
  echo "kept the old default as $KEEP"
elif [[ -L "$DEFAULT" ]]; then
  echo "old default pointed at $(readlink -f "$DEFAULT"); that file is kept"
fi
ln -sfn "$NEW_REAL" "$DEFAULT"
ln -sfn "$NEW_REAL.provenance.json" "$DEFAULT.provenance.json"
echo "default: $DEFAULT -> $NEW_REAL"

if (( VERIFY )); then
  bash "$HERE/../../hpc-hail-job/templates/submit-hail-job.sh" --sif "$DEFAULT" --name "$NAME" \
    --cpus 4 --mem-gb 16 --time 0-00:20:00 --wait -- hvantk utils check-install \
    && echo "PROMOTED IMAGE OK: $DEFAULT" \
    || { echo "ERROR: the promoted image failed the install check; point the default back by hand" >&2; exit 1; }
fi
