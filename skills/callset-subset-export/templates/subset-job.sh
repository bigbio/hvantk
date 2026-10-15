#!/usr/bin/env bash
# Subset a callset and deliver a verified VCF as a chain of four jobs, or list a
# callset's samples as one job. Login-node orchestrator.
#
#   subset-job.sh --list-samples <callset.vds>
#   subset-job.sh --vds <callset.vds> --samples <list> --out <dir> --prefix <name>
#                 [--requested-samples <original list>] [--min-ac 1] [--n-partitions N]
#                 [--cpus 8] [--mem-gb 32] [--time 0-00:30:00] [--partition <P>] [--name NAME]
#
# Each job waits for the previous one; the chain stops at the first failure. Prints the
# job ids, one per line, and leaves <out>/recipe.txt with the exact commands.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SUBMIT="$HERE/../../hpc-hail-job/templates/submit-hail-job.sh"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"; HERE_REL="${HERE#"$R"/}"; T="$HVANTK_REPO/$HERE_REL"

LIST=""; VDS=""; SAMPLES=""; REQUESTED=""; OUT=""; PREFIX="subset"; MIN_AC=1; NPART=""
CPUS=8; MEM_GB=32; TIME="0-00:30:00"; PARTITION=""; NAME=""
while (( $# )); do
  case "$1" in
    --list-samples) LIST="$2"; shift 2 ;;
    --vds) VDS="$2"; shift 2 ;;
    --samples) SAMPLES="$2"; shift 2 ;;
    --requested-samples) REQUESTED="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --min-ac) MIN_AC="$2"; shift 2 ;;
    --n-partitions) NPART="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem-gb) MEM_GB="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    *) echo "ERROR: unknown option $1" >&2; exit 2 ;;
  esac
done
COMMON=(--cpus "$CPUS" --mem-gb "$MEM_GB" --time "$TIME" --wait)
[[ -n "$PARTITION" ]] && COMMON+=(--partition "$PARTITION")
NAME="${NAME:-subset-$PREFIX}"

if [[ -n "$LIST" ]]; then
  [[ "$LIST" == /* ]] || { echo "ERROR: --list-samples needs an absolute path" >&2; exit 2; }
  exec bash "$SUBMIT" "${COMMON[@]}" --name "$NAME-list" -- python "$T/subset-vds.py" --list-samples "$LIST"
fi
[[ -n "$VDS" && -n "$SAMPLES" && -n "$OUT" ]] || { echo "ERROR: --vds, --samples and --out are required" >&2; exit 2; }
[[ "$VDS" == /* && "$SAMPLES" == /* && "$OUT" == /* && ( -z "$REQUESTED" || "$REQUESTED" == /* ) ]] || { echo "ERROR: give absolute paths; jobs run from the log directory" >&2; exit 2; }
mkdir -p "$OUT"
COMMIT="$(git -C "$HVANTK_REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
IMAGE_SHA="$(sha256sum "$(readlink -f "$HVANTK_SIF")" | cut -d' ' -f1)"
NP=(); [[ -n "$NPART" ]] && NP=(--n-partitions "$NPART")
RECIPE="$OUT/recipe.txt"
{
  echo "python subset-vds.py --vds $VDS --samples $SAMPLES --out $OUT/$PREFIX.vds"
  echo "hvantk hgc vds2mt -i $OUT/$PREFIX.vds -o $OUT/$PREFIX.mt ${NP[*]:-}"
  echo "hvantk hgc mt2vcf -i $OUT/$PREFIX.mt -o $OUT/$PREFIX.vcf.bgz --min-ac $MIN_AC"
  echo "python verify-export.py $OUT --prefix $PREFIX --samples $SAMPLES --min-ac $MIN_AC"
} > "$RECIPE"

bash "$SUBMIT" "${COMMON[@]}" --name "$NAME-1subset" -- python "$T/subset-vds.py" --vds "$VDS" --samples "$SAMPLES" --out "$OUT/$PREFIX.vds"
bash "$SUBMIT" "${COMMON[@]}" --name "$NAME-2vds2mt" -- hvantk hgc vds2mt -i "$OUT/$PREFIX.vds" -o "$OUT/$PREFIX.mt" ${NP[@]+"${NP[@]}"}
bash "$SUBMIT" "${COMMON[@]}" --name "$NAME-3mt2vcf" -- hvantk hgc mt2vcf -i "$OUT/$PREFIX.mt" -o "$OUT/$PREFIX.vcf.bgz" --min-ac "$MIN_AC"
REQ=(); [[ -n "$REQUESTED" ]] && REQ=(--requested-samples "$REQUESTED")
# The verify job reads and hashes the whole VCF single-threaded: it gets the chain's time.
VERIFY=(--cpus 2 --mem-gb 8 --time "$TIME" --wait); [[ -n "$PARTITION" ]] && VERIFY+=(--partition "$PARTITION")
bash "$SUBMIT" "${VERIFY[@]}" --name "$NAME-4verify" -- python "$T/verify-export.py" "$OUT" \
  --prefix "$PREFIX" --samples "$SAMPLES" --source "$VDS" --commit "$COMMIT" --image-sha256 "$IMAGE_SHA" --recipe "$RECIPE" --min-ac "$MIN_AC" ${REQ[@]+"${REQ[@]}"}
