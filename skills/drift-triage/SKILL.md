---
name: drift-triage
description: 'Use when hvantk drift or a drift-bot pull request reports a dataset as drifted or its probe as failed: classifying the signal (routine content change, schema change, broken probe), choosing between accepting the new fingerprint and repairing the plugin, rebuilding and verifying the dataset against its committed snapshots, and recording the fingerprint, the rebuild-ledger row and the rebuild. Also when a build that used to work fails with a missing column, a key field not found, or a header marker not found.'
user-invocable: true
argument-hint: "<provider:dataset> [--report drift-report.json] [--raw-dir DIR]"
metadata:
  stage: annotation
  status: provisional
  scale: node
---

# Triage a drift signal and repair the plugin

Read the `hvantk:workflow-conventions` skill (`skills/workflow-conventions/SKILL.md`)
first; this skill assumes its site-profile rule, its section contract and its
measurement contract.

## 1. Scope and outcome

Turn one dataset's drift signal into a reviewable change: the plugin repaired or
left alone with a reason, the dataset rebuilt from the current upstream file and
verified against its own schema snapshot, the fingerprint regenerated, and the
rebuild-ledger row recorded and marked rebuilt.

Use it on a `drifted` report, `probe_failed`, an open drift-bot pull request, or a
build that used to work and now fails with a missing column, a key field not found, or
a header marker not found. Classifying and repairing (steps 2, 4) touch only the
checkout, no image, no job; reproducing, rebuilding, verifying and the bookkeeping
(1, 3, 5, 6) all run through the image, where the fingerprint and ledger files live.

Not this skill: writing a new plugin (a resource spec under `hvantk/skills/<provider>/`),
changing a dataset's output contract on purpose (`references/repair-patterns.md` says
where that line is), building many datasets at once (`dataset-build-at-scale`), or
redesigning a probe beyond making it answer again.

## 2. Preflight

1. Site profile: `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`; stop and ask for any missing value. Needed by every image-touching step (1, 3, 5, 6); classifying (2) and repairing (4) are checkout-only.
2. Checkout: the commit under review (the drift bot's branch, or `dev`) is checked out at `$HVANTK_REPO`, clean tree -- `git -C "$HVANTK_REPO" status --porcelain` -- so the repair you make is the only diff.
3. Image: confirm the path and record its digest (`sha256sum "$(readlink -f "$HVANTK_SIF")"`); rebuild and verify jobs run inside it.
4. Egress: a login node with egress answers `hvantk drift` directly. Without it -- a compute node, a sandbox, or to pin the signal to a file in hand -- `templates/probe-offline.py` answers the same probe locally (`references/bookkeeping.md` says which).
5. Inputs and outputs: a raw directory for the current upstream file, a writable output path for the rebuild.

## 3. Procedure

1. Reproduce the signal, through the image so the checkout's own fingerprint and ledger
   are what gets read (never a bare `hvantk` or `python -m hvantk` here -- the PATH/
   stale-copy trap below):
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
   moving is routine; `probe_error` means the probe is broken, not the dataset --
   three different causes there, each with its own fix (same ref).
3. Get the current upstream file into a dated raw directory and compare its header
   with the committed fingerprint, not the plugin's fixture (GenCC's has 14 of 31
   columns) -- through the image again (downloads, so it also needs egress):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk reprocess <provider:dataset> --raw-dir <dir> --output <dir>/unused.ht \
       --skip-build --no-check-drift --plugin-arg overwrite=true
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     python "$HVANTK_REPO/skills/drift-triage/templates/probe-offline.py" diff \
       --dataset <provider:dataset> --file <new-file> --against-committed
   ```
   Then read the plugin's rename map (`hvantk/skills/<provider>/shared/constants.py`, or
   the builder) to see which declared source columns the new file lacks.
4. Repair per `references/repair-patterns.md`: a renamed column gets its new source
   name added to the rename map, keeping the old one; a removed column feeding an
   output field needs a decision with the user (the output contract changes); an
   added column needs nothing unless it should be mapped; a moved header marker is a
   builder fix.
5. Rebuild as a job, then verify and run the plugin's own tests. The build is the one
   step worth sizing for: submit it in the background and poll, as
   `skills/hpc-hail-job/SKILL.md` section 3 steps 3-4 describe (a foreground `--wait`
   outlives most tool timeouts):
   ```bash
   source "$HVANTK_SITE_PROFILE" && nohup bash -c 'source "$HVANTK_SITE_PROFILE" && \
     bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" --name <run>-build \
       --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
       hvantk reprocess <provider:dataset> --raw-dir <dir> --output <out> \
         --skip-download --no-check-drift' > <run>-build.submit.log 2>&1 &
   ```
   The job id is printed at the top of that log once submitted; poll with `squeue -j
   <id> ...`/`sacct -j <id> ...`. Verifying and testing are quick enough to `--wait` on:
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
   `verify-builds.py` is `dataset-build-at-scale`'s template, referenced here by path.
6. Bookkeeping, per `references/bookkeeping.md` -- regenerate BEFORE marking rebuilt
   (the reverse order is a trap, below):
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
   `resources/drift_ledger.json`, and the fixture or snapshots only if intentionally
   regenerated (say so in the commit). Commit as `fix(<provider>): accept <what
   changed upstream>`.

## 4. Scale-up

- A routine batch (several datasets drifted together): accept every one, rebuild as
  one array through `dataset-build-at-scale`, mark each rebuilt individually.
- A large, variant-level dataset: accept the fingerprint now, schedule the rebuild
  separately; `rebuilt_at: null` keeps the debt visible until it runs.
- A slow or flaky probe: `hvantk drift --timeout <seconds>`.
- A schema signal shared by several datasets on one baseline (one probe, one signal,
  per `hvantk/skills/_conventions/SKILL.md` section 12): repair once, then rebuild
  every dataset sharing it, not only the one named in the report.

## 5. Traps

- `probe_failed`, by cause (`probe_error` tells them apart; full table in
  `references/triage.md`). Network error: a login node with egress, or answer
  `templates/probe-offline.py` from a downloaded copy. No baseline yet, or one shaped
  like a placeholder (empty checksum, sentinel string, epoch `fetched_at`):
  `--regenerate`. The probe's own code raising on content it DID fetch (a missing
  marker, say): repair `drift_probe.py`; bump `PROBE_VERSION` only for an
  OUTPUT-shape change (`references/bookkeeping.md`).
- A merged drift-bot PR, but the build is still broken. Cause: the bot only bumps the
  fingerprint and ledger row, never a builder. Fix: run this skill's repair step
  regardless of which PR merged it.
- A hand-edited `drift_fingerprint.json`. Cause: only a placeholder-shaped baseline
  (trap above) is rejected outright; any other edit is accepted as a file but then
  never matches a live probe, so the dataset reports `drifted` forever. Fix:
  `--regenerate`, never by hand.
- The printed `observed`/`expected` differ in `fetched_at` or `informational`, which
  looks like a reason to act. Cause: both show in full but are excluded from `diff`
  before it is even computed. Fix: read `diff.changed`, never the raw blocks.
- `--regenerate` run minutes after the report shows a different `content_length`.
  Cause: upstream moved again in the gap; the bot skips this case rather than commit
  an unvalidated result. Fix: use the offline form pinned to the file you built from.
- No ledger row appears after `--regenerate`. Cause: a first capture, or a
  `probe_version`-only change -- neither is an upstream change by itself. Fix: add
  the row by hand if upstream genuinely moved too.
- `--mark-rebuilt` run BEFORE `--regenerate`. Cause: `--regenerate` REPLACES the row
  outright (`rebuilt_at: null`), so running it first is just overwritten -- order
  alone is not dangerous. The real false-clean case: `--regenerate` records NO row (a
  first capture or `probe_version`-only change) while an older row already exists, so
  `--mark-rebuilt` marks THAT one instead and `--ledger` reports clean. Fix:
  regenerate first, always; confirm `last_upstream_change` is from THIS run.
- An added upstream column missing where expected, or an incomplete-looking schema
  snapshot. Cause: a rename map only renames columns it recognizes; anything else
  passes through under its upstream name, so a snapshot (a small fixture) is a FLOOR,
  not an exact description. Fix: compare the contract (key, declared fields), not
  equality; compare values against your own build, never the snapshot.
- Editing a snapshot by hand to make a build pass. Cause: this hides a real contract
  change. Fix: regenerate snapshots only for an intended change
  (`references/repair-patterns.md`), and say so in the commit.
- Plugin tests fail with a fixture-not-found error. Cause: fixture paths are relative
  to the checkout root. Fix: `templates/run-plugin-tests.py` chdirs there first --
  run it, not pytest.
- A bare `hvantk`/`python -m hvantk`, or a bookkeeping syntax error. Cause: outside
  the image `hvantk` is not on PATH and the login node's Python predates these
  templates; `-m hvantk`/`-c` also put the CWD first on `sys.path`, silently reading
  the image's stale copy without the checkout on PYTHONPATH. Fix: route every
  command through `run-in-image.sh`, as every step above does.

## 6. Verification

1. The receipt from `dataset-build-at-scale/templates/verify-builds.py`: the schema
   satisfies the snapshot as a floor (every declared field present with its type;
   extra upstream fields expected) and the row count as expected -- compared against
   your own previous build, never the fixture's hand-authored keys.
2. Through the image as in step 1 (or the offline form against the file just
   validated): `hvantk drift <provider:dataset>` reports clean, and `hvantk drift
   --ledger` no longer lists the dataset.
3. The dataset's ledger row: `signal`, a non-empty `accepted_in`, and
   `last_upstream_change`/`rebuilt_at` from THIS run -- not a stale row the
   `--mark-rebuilt`-before-`--regenerate` trap left behind.
4. The plugin's own tests green: `run-plugin-tests.py --repo "$HVANTK_REPO" --dataset
   <provider:dataset>`, run as a job like step 5's verify.
5. The run notes next to the rebuild: commit, container digest, and the job id of
   every job submitted.
