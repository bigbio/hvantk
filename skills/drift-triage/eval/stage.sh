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

# Exactly one .csv and nothing else in a raw directory: a second, stale file (last
# run's download left behind, say) is precisely the "found 2" raw_inputs.py trap, and
# it needs to be caught here rather than confusing a build later.
assert_exactly_one_csv() {
  local dir="$1" n_csv n_total
  n_csv="$(find "$dir" -maxdepth 1 -name '*.csv' -type f | wc -l | tr -d ' ')"
  n_total="$(find "$dir" -maxdepth 1 -type f | wc -l | tr -d ' ')"
  [[ "$n_csv" == "1" && "$n_total" == "1" ]] \
    || { echo "ERROR: expected exactly one .csv file and nothing else under $dir; found $n_total file(s), $n_csv matching *.csv" >&2; exit 1; }
}

EXISTING="$(find "$ORIGINAL_DIR" -maxdepth 1 -name '*.csv' -print -quit)"
if [[ -n "$EXISTING" ]]; then
  echo "already staged: $EXISTING"
else
  # The plugin's own downloader from the checkout, so the URL has one home
  # (hvantk/skills/clingen/shared) and the checkout's code is what runs.
  "$CONTAINER_BIN" exec "${BINDS[@]}" --env PYTHONPATH="$HVANTK_REPO" "$HVANTK_SIF" \
    hvantk download clingen --output-dir "$ORIGINAL_DIR"
fi
assert_exactly_one_csv "$ORIGINAL_DIR"
ORIGINAL="$(find "$ORIGINAL_DIR" -maxdepth 1 -name '*.csv' -print -quit)"
MUTATED="$MUTATED_DIR/$(basename "$ORIGINAL")"

# Simulate the upstream rename this eval is built around: the key column's header,
# "GENE ID (HGNC)", becomes "HGNC ID". Only the header line (the one starting with
# "GENE SYMBOL") is touched, and only its first match -- sed's s/// with no g flag.
sed '/^"GENE SYMBOL"/ s/"GENE ID (HGNC)"/"HGNC ID"/' "$ORIGINAL" > "$MUTATED"
assert_exactly_one_csv "$MUTATED_DIR"

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

# Python-3.6-safe: this runs with the login node's own python3. Also relativizes
# fingerprint_path (an absolute path under this cluster's checkout) and probe_ref
# (never a path today, handled the same way for uniformity) so the staged report --
# which the task hands to an agent -- never names where the checkout happens to live.
python3 - "$REPORT" "$HVANTK_REPO" <<'PY'
import json
import sys

report_path, repo = sys.argv[1], sys.argv[2]
prefix = repo.rstrip("/") + "/"

report = json.load(open(report_path))
rows = report if isinstance(report, list) else [report]
for row in rows:
    for field in ("fingerprint_path", "probe_ref"):
        value = row.get(field)
        if isinstance(value, str) and value.startswith(prefix):
            row[field] = value[len(prefix):]

row = rows[0]
if row.get("status") != "drifted":
    sys.exit("expected status drifted, got %r" % (row.get("status"),))
changed = (row.get("diff") or {}).get("changed") or {}
if "headers" not in changed:
    sys.exit("expected 'headers' among the changed keys, got %r" % (sorted(changed),))

with open(report_path, "w") as fh:
    fh.write(json.dumps(report, indent=2) + "\n")
print("drift report OK: status=drifted, headers changed; paths relativized to the checkout")
PY

( cd "$STAGE" && sha256sum raw-original/clingen/* raw/clingen/* drift-report.json > staged.sha256 )
cat "$STAGE/staged.sha256"
echo "staged under $STAGE"
