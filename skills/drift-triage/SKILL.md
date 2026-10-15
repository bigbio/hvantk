---
name: drift-triage
description: 'Use when hvantk drift or a drift-bot pull request reports a dataset as drifted or its probe as failed: classifying the signal (routine content change, schema change, broken probe), choosing between accepting the new fingerprint and repairing the plugin, rebuilding and verifying the dataset against its committed snapshots, and recording the fingerprint, the rebuild-ledger row and the rebuild. Also when a build that used to work fails with a missing column, a key field not found, or a header marker not found.'
user-invocable: true
argument-hint: "<provider:dataset> [--report drift-report.json] [--raw-dir DIR]"
metadata:
  stage: annotation
  status: provisional
  scale: laptop
---

# Triage a drift signal and repair the plugin

Read `skills/_conventions/SKILL.md` first; this skill assumes its site-profile rule,
its section contract and its measurement contract.

## 1. Scope and outcome

Turn one dataset's drift signal into a reviewable change: the plugin repaired or
deliberately left alone with a reason, the dataset rebuilt from the current upstream
file and verified against its own committed snapshots, the fingerprint regenerated, and
the rebuild-ledger row recorded and marked rebuilt.

Use it on a `drifted` report, on `probe_failed`, on an open drift-bot pull request, or on
a build that used to work and now fails with a missing column, a key field not found, or
a header marker not found. Classifying and repairing (steps 1-4 below) need nothing but
the checkout and, for step 1's live form, a login node with egress -- no cluster job.
Rebuilding and verifying at hvantk's usual scale (step 5 on) go through the same
container-image jobs every other hvantk workflow skill uses.

Not this skill: writing a new plugin (a resource spec under `hvantk/skills/<provider>/`),
changing a dataset's output contract on purpose (`references/repair-patterns.md` says
where that line is), building many datasets at once (`dataset-build-at-scale`), or
redesigning a probe beyond making it answer again.

## 2. Preflight

1. Site profile: `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`; stop and ask for any missing value. Needed from step 5 on; steps 1-4 touch no profile value.
2. Checkout: the commit under review (the drift bot's branch, or `dev`) is what is checked out at `$HVANTK_REPO`, with a clean tree -- `git -C "$HVANTK_REPO" status --porcelain` -- so the repair you make is the only diff.
3. Image: confirm the image path and record its digest (`sha256sum "$(readlink -f "$HVANTK_SIF")"`); the rebuild and verify jobs run inside it.
4. Egress: a login node with egress answers `hvantk drift` against upstream directly. Without egress -- a compute node, a sandbox -- or to pin the signal to a file already in hand, `templates/probe-offline.py` answers the same probe from a local copy instead; `references/bookkeeping.md` says which situation calls for which.
5. Inputs and outputs: a raw directory for the current upstream file, and a writable output path for the rebuild.

## 3. Procedure

1. Reproduce the signal:
   ```bash
   source "$HVANTK_SITE_PROFILE" && "${HVANTK_CONTAINER_BIN:-singularity}" exec \
     -B "$HVANTK_WORK:$HVANTK_WORK" -B "$HVANTK_REPO:$HVANTK_REPO" \
     --env PYTHONPATH="$HVANTK_REPO" "$HVANTK_SIF" \
     hvantk drift <provider:dataset> --json > drift-report.json
   ```
   With no egress, or to pin the signal to a file already on disk:
   `python "$HVANTK_REPO/skills/drift-triage/templates/probe-offline.py" drift --dataset <provider:dataset> --file <upstream-file> -- --json`.
   A drift-bot pull request's body carries the same diff.
2. Classify with `references/triage.md`: `headers`/`checksums` in `diff.changed`, or a
   top-level key added or removed, is a schema signal; only `extras`/`source_version`
   moving is routine; a `probe_error` means the probe is broken, not that the dataset
   drifted.
3. Get the current upstream file into a dated raw directory and compare its header with
   the plugin's fixture:
   ```bash
   hvantk reprocess <provider:dataset> --raw-dir <dir> --output <dir>/unused.ht \
     --skip-build --no-check-drift --plugin-arg overwrite=true
   python "$HVANTK_REPO/skills/drift-triage/templates/probe-offline.py" diff \
     --dataset <provider:dataset> --file <new-file> --against <fixture>
   ```
   Then read the plugin's rename map (`hvantk/skills/<provider>/shared/constants.py`, or
   the builder) to see which declared source columns the new file lacks.
4. Repair per `references/repair-patterns.md`: a renamed column gets its new source name
   added to the rename map, keeping the old one; a removed column that feeds an output
   field needs a decision with the user, since the output contract changes; an added
   column needs nothing unless wanted; a moved header marker is a builder fix.
5. Rebuild, verify against the committed snapshots, and run the plugin's own tests, each
   as a job:
   ```bash
   bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" --name <run>-build \
     --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     hvantk reprocess <provider:dataset> --raw-dir <dir> --output <out> \
       --skip-download --no-check-drift
   bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" --name <run>-verify \
     --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/dataset-build-at-scale/templates/verify-builds.py" \
       --build <provider:dataset>=<out> --receipt <out>.receipt.json
   bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" --name <run>-test \
     --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/drift-triage/templates/run-plugin-tests.py" \
       --repo "$HVANTK_REPO" --provider <provider>
   ```
   `verify-builds.py` is `dataset-build-at-scale`'s template (`--build DS=OUT --receipt
   PATH`, the same contract its own skill documents); referenced here by path, not
   restated.
6. Bookkeeping, per `references/bookkeeping.md`:
   ```bash
   hvantk drift --regenerate <provider:dataset> --accepted-in "<PR or branch>"
   hvantk drift --mark-rebuilt <provider:dataset>
   hvantk drift --ledger   # must no longer list <provider:dataset>
   ```
   On a machine with no egress, or to pin the baseline to the file just validated, the
   first command becomes `python "$HVANTK_REPO/skills/drift-triage/templates/probe-offline.py"
   drift --dataset <provider:dataset> --file <new-file> -- --regenerate --accepted-in "<PR or branch>"`.
7. Review the diff before committing: the repair, `tests/drift_fingerprint.json`,
   `resources/drift_ledger.json`, and the fixture or snapshots ONLY if regenerated on
   purpose (say so in the commit message). Commit as
   `fix(<provider>): accept <what changed upstream>`.

## 4. Scale-up

- A routine batch (several datasets drifted together): accept every one, rebuild them as
  one array through `dataset-build-at-scale`, then mark each rebuilt individually --
  `--mark-rebuilt` takes one dataset at a time.
- A large, variant-level dataset: accept the fingerprint now and schedule the rebuild
  separately; the ledger's `rebuilt_at: null` keeps the debt visible until it runs.
- A slow or flaky probe: `hvantk drift --timeout <seconds>`.
- A schema signal shared by several datasets on one baseline (one probe, one signal, per
  `hvantk/skills/_conventions/SKILL.md` section 12): repair once, then rebuild every
  dataset that shares it, not only the one named in the report.

## 5. Traps

- `probe_failed` triaged as if it were drift. Cause: the probe itself raised (network,
  parsing, a missing baseline), which carries no information about the dataset. Fix:
  repair `drift_probe.py` and bump `PROBE_VERSION`; no ledger row for that alone.
- A diff that only shows `fetched_at`/`informational` moving. Cause: those keys are
  already excluded from the comparison that produced the report. Fix: there is nothing
  to accept; re-run the check if the report still looks wrong.
- A drift-bot pull request merged, but the build is still broken. Cause: the bot only
  bumps the fingerprint and records the ledger row; it has never repaired a builder. Fix:
  run this skill's repair step regardless of which PR merged the fingerprint.
- A hand-edited `drift_fingerprint.json`. Cause: a baseline not captured from a live
  probe is rejected outright, or never matches a real observation again. Fix:
  `--regenerate`, never a manual edit.
- `--regenerate` run minutes after the report shows a different `content_length`. Cause:
  upstream moved again in the gap between the two calls; the drift bot skips exactly
  this case rather than commit a result nobody validated. Fix: use the offline form
  pinned to the file you actually built and tested from.
- No ledger row appears after `--regenerate`. Cause: a first capture, or a
  `probe_version`-only change -- neither is an upstream change by itself. Fix: add the
  row by hand if upstream genuinely also moved.
- `--mark-rebuilt` refuses the dataset. Cause: it has no ledger row yet. Fix:
  `--regenerate` first; `--mark-rebuilt` only clears a row that exists.
- A rename-map builder stays quiet on a column it no longer recognizes. Cause: some
  builders only warn and let the output schema silently gain the new name. Fix: check
  the schema snapshot after a rebuild, not just the exit code.
- Editing `sample_rows.json` by hand to make a build pass. Cause: this hides a real
  contract change instead of fixing it. Fix: `--regenerate-snapshots` only for an
  intended change, and say so in the commit.
- Plugin tests fail with a fixture-not-found error. Cause: fixture paths are relative to
  the checkout root. Fix: `templates/run-plugin-tests.py` chdirs there first; run it
  rather than invoking pytest from elsewhere.
- A bookkeeping command raises a syntax error on the login node. Cause: its Python
  predates what these scripts are written against. Fix: run every check through the
  image, as step 5 and 6 already do.

## 6. Verification

1. The receipt `dataset-build-at-scale/templates/verify-builds.py` wrote: schema match,
   row count, every sample key present.
2. `hvantk drift <provider:dataset>` (or the offline form against the file just
   validated) reports clean.
3. `hvantk drift --ledger` no longer lists the dataset.
4. The dataset's ledger row: `signal`, a non-empty `accepted_in`, and `rebuilt_at` no
   earlier than `last_upstream_change`.
5. The plugin's own tests green: `templates/run-plugin-tests.py --repo "$HVANTK_REPO"
   --provider <provider>`.
6. The run notes next to the rebuild: commit, container digest, and the job id of every
   job submitted.
