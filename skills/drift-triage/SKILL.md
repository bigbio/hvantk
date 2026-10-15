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
file and verified against its own committed schema snapshot, the fingerprint
regenerated, and the rebuild-ledger row recorded and marked rebuilt.

Use it on a `drifted` report, on `probe_failed`, on an open drift-bot pull request, or on
a build that used to work and now fails with a missing column, a key field not found, or
a header marker not found. Classifying and repairing (steps 2 and 4 below) touch only the
checkout, no image and no job; reproducing the signal, rebuilding, verifying and the
bookkeeping (1, 3, 5, 6) all run through the container image, because that is where the
fingerprint and ledger files that get written actually live.

Not this skill: writing a new plugin (a resource spec under `hvantk/skills/<provider>/`),
changing a dataset's output contract on purpose (`references/repair-patterns.md` says
where that line is), building many datasets at once (`dataset-build-at-scale`), or
redesigning a probe beyond making it answer again.

## 2. Preflight

1. Site profile: `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`; stop and ask for any missing value. Needed by every step that touches the image (1, 3, 5, 6); classifying (2) and repairing (4) are checkout-only work.
2. Checkout: the commit under review (the drift bot's branch, or `dev`) is what is checked out at `$HVANTK_REPO`, with a clean tree -- `git -C "$HVANTK_REPO" status --porcelain` -- so the repair you make is the only diff.
3. Image: confirm the image path and record its digest (`sha256sum "$(readlink -f "$HVANTK_SIF")"`); the rebuild and verify jobs run inside it.
4. Egress: a login node with egress answers `hvantk drift` against upstream directly. Without egress -- a compute node, a sandbox -- or to pin the signal to a file already in hand, `templates/probe-offline.py` answers the same probe from a local copy instead; `references/bookkeeping.md` says which situation calls for which.
5. Inputs and outputs: a raw directory for the current upstream file, and a writable output path for the rebuild.

## 3. Procedure

1. Reproduce the signal, through the image so the checkout's own fingerprint and ledger
   are what gets read (never a bare `hvantk` or `python -m hvantk` here: outside the
   image it is either not on PATH at all, or its Python is too old for these templates;
   inside the image without the checkout on PYTHONPATH it reads the image's own, stale
   copy and cannot write back to this one):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk drift <provider:dataset> --json > drift-report.json
   ```
   With no egress, or to pin the signal to a file already on disk:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     python "$HVANTK_REPO/skills/drift-triage/templates/probe-offline.py" drift \
       --dataset <provider:dataset> --file <upstream-file> -- --json
   ```
   A drift-bot pull request's body carries the same diff.
2. Classify with `references/triage.md`: `headers`/`checksums` in `diff.changed`, or a
   top-level key added or removed, is a schema signal; only `extras`/`source_version`
   moving is routine; a `probe_error` means the probe is broken, not that the dataset
   drifted -- and `references/triage.md` says which of three different causes that is,
   because the fix differs per cause.
3. Get the current upstream file into a dated raw directory and compare its header with
   the plugin's fixture, through the image again (this one downloads, so it also needs
   egress):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk reprocess <provider:dataset> --raw-dir <dir> --output <dir>/unused.ht \
       --skip-build --no-check-drift --plugin-arg overwrite=true
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     python "$HVANTK_REPO/skills/drift-triage/templates/probe-offline.py" diff \
       --dataset <provider:dataset> --file <new-file> --against <fixture>
   ```
   Then read the plugin's rename map (`hvantk/skills/<provider>/shared/constants.py`, or
   the builder) to see which declared source columns the new file lacks.
4. Repair per `references/repair-patterns.md`: a renamed column gets its new source name
   added to the rename map, keeping the old one; a removed column that feeds an output
   field needs a decision with the user, since the output contract changes; an added
   column needs nothing unless it should be mapped; a moved header marker is a builder
   fix.
5. Rebuild as a job, then verify and run the plugin's own tests. The build is the one
   step worth sizing for: submit it in the background and poll, exactly as
   `skills/hpc-hail-job/SKILL.md` section 3 steps 3-4 describe (a foreground `--wait`
   outlives most tool timeouts):
   ```bash
   source "$HVANTK_SITE_PROFILE" && nohup bash -c 'source "$HVANTK_SITE_PROFILE" && \
     bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" --name <run>-build \
       --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
       hvantk reprocess <provider:dataset> --raw-dir <dir> --output <out> \
         --skip-download --no-check-drift' > <run>-build.submit.log 2>&1 &
   ```
   The job id is printed at the top of that log once submitted; poll it with `squeue -j
   <id> ...` / `sacct -j <id> ...` as the same section shows. Verifying and testing are
   quick enough for a foreground `--wait`:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" --name <run>-verify \
     --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/dataset-build-at-scale/templates/verify-builds.py" \
       --build <provider:dataset>=<out> --receipt <out>.receipt.json
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" --name <run>-test \
     --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/drift-triage/templates/run-plugin-tests.py" \
       --repo "$HVANTK_REPO" --dataset <provider:dataset>
   ```
   `verify-builds.py` is `dataset-build-at-scale`'s template (`--build DS=OUT --receipt
   PATH`); referenced here by path, not restated.
6. Bookkeeping, per `references/bookkeeping.md` -- regenerate BEFORE marking rebuilt; the
   reverse order is a trap (below), not just a style preference:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk drift --regenerate <provider:dataset> --accepted-in "<PR or branch>"
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk drift --mark-rebuilt <provider:dataset>
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk drift --ledger   # must no longer list <provider:dataset>
   ```
   On a machine with no egress, or to pin the baseline to the file just validated, the
   first command becomes:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     python "$HVANTK_REPO/skills/drift-triage/templates/probe-offline.py" drift \
       --dataset <provider:dataset> --file <new-file> -- --regenerate --accepted-in "<PR or branch>"
   ```
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

- `probe_failed` from a network error. Cause: no egress from here. Fix: a login node
  with egress, or `templates/probe-offline.py` against a downloaded copy.
- `probe_failed` from a missing or hand-seeded baseline. Cause: no committed fingerprint
  yet, or one shaped like a stub rather than a live capture. Fix: `--regenerate`.
- `probe_failed` from the probe itself raising. Cause: its own code is broken (a parse
  error, a changed response shape) -- not the dataset. Fix: repair `drift_probe.py`;
  bump `PROBE_VERSION` only when its OUTPUT shape changed, knowing that ANY version
  difference suppresses `--regenerate`'s ledger row (`references/bookkeeping.md`), so
  bumping it for something unrelated quietly loses the next real signal's row.
- A drift-bot pull request merged, but the build is still broken. Cause: the bot only
  bumps the fingerprint and records the ledger row; it has never repaired a builder. Fix:
  run this skill's repair step regardless of which PR merged the fingerprint.
- A hand-edited `drift_fingerprint.json`. Cause: only a baseline shaped like a stub or a
  never-captured placeholder is rejected outright; any other hand-edit is accepted as a
  file but then simply never matches a live probe again, so the dataset reports
  `drifted` forever. Fix: `--regenerate`, never a manual edit.
- The printed `observed`/`expected` differ in `fetched_at` or under `informational`, and
  that looks like a reason to act. Cause: both are shown in full, but excluded from
  `diff` before it is even computed -- `informational` never appears in a diff at all.
  Fix: read `diff.changed`, not the raw `observed`/`expected` blocks.
- `--regenerate` run minutes after the report shows a different `content_length`. Cause:
  upstream moved again in the gap between the two calls; the drift bot skips exactly
  this case rather than commit a result nobody validated. Fix: use the offline form
  pinned to the file you actually built and tested from.
- No ledger row appears after `--regenerate`. Cause: a first capture, or a
  `probe_version`-only change -- neither is an upstream change by itself. Fix: add the
  row by hand if upstream genuinely also moved in the same window.
- `--mark-rebuilt` run BEFORE `--regenerate` marks the wrong row. Cause: if the dataset
  already has an older row from an earlier, unrelated signal, `--mark-rebuilt` happily
  marks THAT one rebuilt, and `--ledger` then reports clean even though THIS signal was
  never recorded -- it only refuses outright when there is no row at all yet. Fix:
  regenerate first, always, then check the row's `last_upstream_change` is this run's,
  not a stale one.
- An added upstream column does not show up where expected, or a schema snapshot looks
  incomplete. Cause: a builder's rename map only renames the columns it recognizes --
  anything else passes through under its own upstream name, so the committed snapshot
  (built from a small fixture) is a FLOOR, not an exact description of a full build. Fix:
  compare against the contract (key present, every declared field present), not for
  exact equality; compare values against your own previous build, not the snapshot.
- Editing a snapshot file by hand to make a build pass. Cause: this hides a real
  contract change instead of fixing it. Fix: regenerate snapshots only for an intended
  change (`references/repair-patterns.md`), and say so in the commit.
- Plugin tests fail with a fixture-not-found error. Cause: fixture paths are relative to
  the checkout root. Fix: `templates/run-plugin-tests.py` chdirs there first; run it
  rather than invoking pytest from elsewhere.
- A bookkeeping command raises a syntax error. Cause: the login node's own Python
  predates what these templates are written against. Fix: every step above already
  runs through the image for exactly this reason; do not shortcut it to save a step.

## 6. Verification

1. The receipt from `dataset-build-at-scale/templates/verify-builds.py`: the schema
   satisfies the committed snapshot as a floor (key matches, every declared field
   present with its declared type -- extra fields from upstream are expected, not a
   failure) and the row count is what you expect. Compare values against your own
   previous build, not against any "sample row" promise -- the fixture's sample keys
   are hand-authored and are not guaranteed to exist in a real build at all.
2. `hvantk drift <provider:dataset>`, through the image as in step 1 (or the offline
   form against the file just validated), reports clean.
3. `hvantk drift --ledger`, through the image the same way, no longer lists the
   dataset.
4. The dataset's ledger row: `signal`, a non-empty `accepted_in`, and
   `last_upstream_change` / `rebuilt_at` that belong to THIS run, not a stale row the
   mark-rebuilt/regenerate ordering trap left behind.
5. The plugin's own tests green: `run-plugin-tests.py --repo "$HVANTK_REPO" --dataset
   <provider:dataset>`, submitted as a job the same way step 5's verify job is.
6. The run notes next to the rebuild: commit, container digest, and the job id of every
   job submitted.
