#!/usr/bin/env bash
# Stage the public input for the drift-triage eval, once per cluster. Login-node work:
# download the live ClinGen file, mutate a copy's header the way this eval's upstream
# change is modelled (the key column renamed), then pin a drift report to that mutated
# copy so the measured task runs with no network. Writes staged.sha256 next to the
# staged files.
#
#   HVANTK_SITE_PROFILE=... bash skills/drift-triage/eval/stage.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/drift-triage"
ORIGINAL_DIR="$STAGE/raw-original/clingen"
MUTATED_DIR="$STAGE/raw/clingen"
mkdir -p "$ORIGINAL_DIR" "$MUTATED_DIR"
CONTAINER_BIN="${HVANTK_CONTAINER_BIN:-singularity}"
BINDS=(-B "$HVANTK_WORK:$HVANTK_WORK")
for p in "$HVANTK_REPO" "$STAGE"; do
  case "$p" in "$HVANTK_WORK"|"$HVANTK_WORK"/*) ;; *) BINDS+=(-B "$p:$p") ;; esac
done
PROBE_OFFLINE="$HVANTK_REPO/skills/drift-triage/templates/probe-offline.py"

EXISTING="$(find "$ORIGINAL_DIR" -maxdepth 1 -name '*.csv' -print -quit)"
if [[ -n "$EXISTING" ]]; then
  echo "already staged: $EXISTING"
else
  # The plugin's own downloader from the checkout, so the URL has one home
  # (hvantk/skills/clingen/shared) and the checkout's code is what runs.
  "$CONTAINER_BIN" exec "${BINDS[@]}" --env PYTHONPATH="$HVANTK_REPO" "$HVANTK_SIF" \
    hvantk download clingen --output-dir "$ORIGINAL_DIR"
fi
ORIGINAL="$(find "$ORIGINAL_DIR" -maxdepth 1 -name '*.csv' -print -quit)"
[[ -n "$ORIGINAL" ]] || { echo "ERROR: no CSV staged under $ORIGINAL_DIR" >&2; exit 1; }
MUTATED="$MUTATED_DIR/$(basename "$ORIGINAL")"

# Simulate the upstream rename this eval is built around: the key column's header,
# "GENE ID (HGNC)", becomes "HGNC ID". Only the header line (the one starting with
# "GENE SYMBOL") is touched, and only its first match -- sed's s/// with no g flag.
sed '/^"GENE SYMBOL"/ s/"GENE ID (HGNC)"/"HGNC ID"/' "$ORIGINAL" > "$MUTATED"

# diff exits 1 when the files differ, which is the expected case here -- capture its
# output rather than letting `set -e`/pipefail treat that as a script failure.
DIFF_OUT="$(diff "$ORIGINAL" "$MUTATED" || true)"
ADDED="$(printf '%s\n' "$DIFF_OUT" | grep -c '^>' || true)"
REMOVED="$(printf '%s\n' "$DIFF_OUT" | grep -c '^<' || true)"
[[ "$ADDED" == "1" && "$REMOVED" == "1" ]] \
  || { echo "ERROR: expected the mutation to change exactly one line (1 removed, 1 added), got removed=$REMOVED added=$ADDED" >&2; exit 1; }

REPORT="$STAGE/drift-report.json"
RC=0
"$CONTAINER_BIN" exec "${BINDS[@]}" --env PYTHONPATH="$HVANTK_REPO" "$HVANTK_SIF" \
  python "$PROBE_OFFLINE" drift --dataset clingen:gene-disease --file "$MUTATED" -- --json \
  > "$REPORT" || RC=$?
[[ "$RC" == "1" ]] \
  || { echo "ERROR: expected the drift probe to report drifted (exit 1), got exit $RC" >&2; exit 1; }

# Python-3.6-safe: this runs with the login node's own python3.
python3 - "$REPORT" <<'PY'
import json
import sys

report = json.load(open(sys.argv[1]))
row = report[0] if isinstance(report, list) else report
if row.get("status") != "drifted":
    sys.exit("expected status drifted, got %r" % (row.get("status"),))
changed = (row.get("diff") or {}).get("changed") or {}
if "headers" not in changed:
    sys.exit("expected 'headers' among the changed keys, got %r" % (sorted(changed),))
print("drift report OK: status=drifted, headers changed")
PY

( cd "$STAGE" && sha256sum raw-original/clingen/* raw/clingen/* drift-report.json > staged.sha256 )
cat "$STAGE/staged.sha256"
echo "staged under $STAGE"
