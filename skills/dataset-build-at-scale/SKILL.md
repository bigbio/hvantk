---
name: dataset-build-at-scale
description: 'Use when building one or many hvantk datasets with hvantk reprocess on a Slurm cluster from staged raw files: raw-directory conventions, staging downloads on the login node, fanning independent builds out as one array job, verifying every output against its plugin''s committed schema and sample-row snapshots, and checking upstream drift before a build. Also when a build fails with "expected exactly one file matching", "has no lifecycle.download declared", an HTTP error on a compute node, or an output-extension error.'
user-invocable: true
argument-hint: "<provider:dataset ...> [--raw-root DIR --out-root DIR --cpus N --mem-gb G --time D-HH:MM:SS]"
metadata:
  stage: annotation
  status: provisional
  scale: node
---

# Build hvantk datasets at scale

Read `skills/_conventions/SKILL.md` first; this skill assumes its site-profile rule,
its section contract and its measurement contract.

## 1. Scope and outcome

Building one or many `hvantk reprocess` datasets from raw files already staged on disk,
as Slurm jobs, with every output verified against its plugin's own committed contract.
Several independent datasets build as one array job rather than as separate
submissions. The outcome, per dataset: the output directory (`.ht`, `.mt` or `.h5ad`)
with `_SUCCESS` and its `.provenance.json` sidecar, a verification receipt, and run
notes carrying the array (or single) job id.

Use it to build any plugin dataset — one, or a batch of several — once its raw input is
on disk or downloadable. Not this skill: a single build with nothing to batch (plain
`hpc-hail-job` is enough on its own); writing or repairing a plugin (the resource specs
under `hvantk/skills/<provider>/SKILL.md`); triaging a drift signal or a broken probe
(`drift-triage`).

## 2. Preflight

1. Site profile, checkout, image, inputs and outputs: the contract's four steps
   (`skills/_conventions/SKILL.md` section 2 / the template in section 7).
2. List what you are building: `hvantk plugins list` and `hvantk plugins describe
   <provider>` through the image say each dataset's `artifact_type` (which decides the
   output extension: `.ht` for a Hail `AnnotationTable`, `.mt` for a `VariantMatrix`,
   `.h5ad` for an `ExpressionMatrix`, `.geneset.json` for a `GeneSet`), whether it
   declares `lifecycle.parse` (so it needs `--intermediate`), and whether
   `acquisition.mode: byo` means you must supply the raw files yourself.
3. Egress: whatever the site profile's `HVANTK_EGRESS` says, downloads happen on the
   login node (`references/raw-dir-conventions.md`); the build jobs themselves always
   run with `--skip-download`, so egress changes who runs the download command and
   when, not the procedure below.

## 3. Procedure

1. **Lay out one raw directory per dataset**: `$HVANTK_WORK/raw/<provider>[/<dataset>]/`,
   holding only that dataset's files. Full rules, including dated files and
   `latest_dated`, in `references/raw-dir-conventions.md`.
2. **Stage the raw files on the login node, through the image**, uniformly for any
   dataset with a downloader:
   ```bash
   source "$HVANTK_SITE_PROFILE" && "${HVANTK_CONTAINER_BIN:-singularity}" exec \
     -B "$HVANTK_WORK:$HVANTK_WORK" -B "$HVANTK_REPO:$HVANTK_REPO" --env PYTHONPATH="$HVANTK_REPO" \
     "$HVANTK_SIF" hvantk reprocess provider:dataset --raw-dir "$HVANTK_WORK/raw/provider" \
       --output "$HVANTK_WORK/raw/provider/unused.ht" --skip-build --no-check-drift
   ```
   add `--plugin-arg overwrite=true` to refresh a stale file. Then check upstream has
   not drifted against the committed fingerprint before spending a job on the build —
   `hvantk drift provider:dataset` through the same `exec` line. A `drifted` result
   whose diff touches `headers` or `checksums` is a schema signal: stop and use
   `drift-triage` first. A clean or `routine` result means it is safe to proceed.
3. **Write the task list**: a TSV, one build per line, absolute paths —
   `<provider:dataset>\t<raw dir>\t<output path>[\t<extra reprocess options>]`. Extra
   options (e.g. `--intermediate /abs/path`, `--plugin-arg reference_genome=GRCh38`) are
   word-split exactly as written.
4. **Submit the array**:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/dataset-build-at-scale/templates/submit-build-array.sh" \
     --list "$HVANTK_WORK/out/tasks.tsv" --name run-build --cpus 4 --mem-gb 16 \
     --time 0-00:30:00 --max-concurrent 3 --dry-run
   ```
   then without `--dry-run`. It prints the array job id; write it into your run notes
   at once as `job: <id>`, next to `commit:` and `container:` from the preflight.
5. **Wait and triage.** Poll `sacct -j <id> -X -o JobID,JobName%30,State,Elapsed`, or add
   `--wait` and run it in the background (a foreground wait outlives most tool
   timeouts). Read each task's log at `<name>_<array id>_<task>.out`. One failed task
   does not fail the array and does not stop the others; resubmit only the failed
   indices with `submit-build-array.sh --list <same list> --tasks i,j ...` under a new
   job name.
6. **Verify every output as a job**:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name run-verify --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/dataset-build-at-scale/templates/verify-builds.py" \
       --list "$HVANTK_WORK/out/tasks.tsv" --receipt "$HVANTK_WORK/out/verify.json"
   ```
   Details in `references/verification.md`.
7. **Finish the run notes.** If a build refreshed a dataset whose rebuild-ledger row
   (`hvantk drift --ledger`) was pending, run `hvantk drift --mark-rebuilt
   provider:dataset` on the login node and commit the ledger alongside the change that
   motivated the rebuild.

## 4. Scale-up

`references/sizing-by-dataset.md` has the full table. In short: one request sizes every
task in the array, so size it for the list's largest dataset — a dataset an
order of magnitude bigger than the rest (dbNSFP next to HGNC, say) belongs in its own
`hpc-hail-job` submission, not in the array. `--max-concurrent` caps how much of the
node share the array can claim at once, independent of how many lines are in the list.
Each task gets its own node-local scratch from `hail-job.sbatch`, so concurrent tasks
never collide even on the same node. Outputs on the parallel filesystem are directories
of many part files; keep each dataset's dated raw file and its `.provenance.json`
sidecar discoverable next to its output. Re-running only the failed indices (`--tasks`)
is the normal way to recover from a partial array failure, under a fresh job id rather
than a retry of the old one.

## 5. Traps

- `expected exactly one file matching ...`. Cause: two files in one raw dir (a stray
  second download, or a parse-stage intermediate written inside it). Fix: one dataset
  per raw directory; give `--intermediate` a path outside the raw dir.
- An HTTP or connection error inside a build job. Cause: compute nodes have no egress.
  Fix: stage the raw file on the login node first; build jobs run with
  `--skip-download`.
- `... has no lifecycle.download declared`. Cause: a BYO dataset (`acquisition.mode:
  byo`). Fix: place the raw file yourself; `hvantk plugins describe <provider>` prints
  the instructions.
- The output-extension error (`<dataset> writes .ht ...; got --output ...`). Cause: the
  artifact type decides the extension (`.ht` table, `.mt` variant matrix, `.h5ad`
  expression), and a mismatch is rejected before the build runs. Fix: match `--output`'s
  extension to what `hvantk plugins describe` names as the artifact type.
- A `--plugin-arg` value arriving as the wrong type. Cause: the coercion rules
  (`true`/`false` to bool, digits to int, a comma to a list) — a single value that
  itself contains a comma cannot be expressed this way. Fix: check `hvantk plugins
  describe` for the parameter's expected type.
- `drift: drifted` printed by the login-node staging step. Cause: that is the probe
  reporting on upstream, not the build failing — `reprocess` never fails on drift. Fix:
  read the diff; a `schema` signal is worth stopping for (`drift-triage`), a `routine`
  one is informational. Build jobs run with `--no-check-drift`: they have no egress to
  probe with anyway.
- `sacct` showing one row `<A>_[2-3]` for a while after submission. Cause: tasks still
  pending can be reported collapsed into one range row. Fix: wait for every task to
  reach a terminal state, not for that row to disappear.
- One failed task does not fail the array. Cause: Slurm array tasks are independent
  jobs sharing a job id. Fix: check per-task states (`sacct -j <A> -X -n -o
  JobID,State`), never the submission's own exit code alone.
- Sample-row values differing from the committed `sample_rows.json` on a full build.
  Cause: the committed snapshot is fixture-derived (dates, array widths can legitimately
  differ at full scale). Fix: `verify-builds.py` checks the schema exactly and the
  sample *keys*' presence; compare row *values* across your own builds, not against the
  fixture snapshot.
- A relative path in the task list. Cause: an array task's working directory is the log
  directory, not the checkout. Fix: absolute paths only, both columns.
- `hvantk: command not found` outside the image. Fix: `python -m hvantk`.
- An output moved without its sidecar. Cause: `.provenance.json` sits next to the
  output, not inside it. Fix: move both together, or neither.
- A task list with a stale raw directory. Cause: yesterday's dated file next to today's
  is tolerated for `latest_dated` plugins (ClinGen, GenCC) but not for one without it
  (HGNC). Fix: keep exactly one matching file per raw directory
  (`references/raw-dir-conventions.md`).

## 6. Verification

1. `sacct -j <A> -X -n -o JobID,State` — every task row `COMPLETED`; `<A>` and `<A>_<i>`
   both count as the array's own ids.
2. `verify-builds.py`'s receipt (`references/verification.md`): per dataset,
   `schema_matches` true, `sample_keys_present` true, the row count and sample rows
   recorded for your own comparison against the next build, and `provenance_sha256`
   tying the verdict to exactly this output.
3. Each output carries `_SUCCESS` and a `.provenance.json` sidecar naming the
   `schema_id` you expected.
4. Run notes carry `commit:`, `container:`, and one `job:` line per job submitted (the
   array counts as one `job:` line; the verification job is another).
