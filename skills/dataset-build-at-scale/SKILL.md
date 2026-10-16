---
name: dataset-build-at-scale
description: 'Use when building one or many hvantk datasets with hvantk reprocess on a Slurm cluster from staged raw files: raw-directory conventions, staging downloads on the login node, fanning independent builds out as one array job, verifying every output against its plugin''s committed schema contract, and checking upstream drift before a build. Also when a build fails with "expected exactly one file matching", "requires data you supply yourself", an HTTP error on a compute node, or an output-extension error.'
user-invocable: true
argument-hint: "<provider:dataset ...> [--list tasks.tsv --cpus N --mem-gb G --time D-HH:MM:SS --max-concurrent K]"
metadata:
  stage: annotation
  status: provisional
  scale: node
---

# Build hvantk datasets at scale

Read the `hvantk:workflow-conventions` skill (`skills/workflow-conventions/SKILL.md`)
first; this skill assumes its site-profile rule, its section contract and its
measurement contract.

## 1. Scope and outcome

Building and verifying `hvantk reprocess` datasets from raw files already staged on
disk, as Slurm jobs. The task-list and array-job machinery is overkill for exactly one
dataset (`hpc-hail-job` alone covers that), but the raw-dir conventions and the
verification contract below are the same whether it is one or several; the array is
how several independent builds share one Slurm request instead of one submission each.
The outcome, per dataset: the output (a `.ht`/`.mt` directory with `_SUCCESS`, or a
single `.h5ad`/`.parquet` file) with its `.provenance.json` sidecar, a verification
receipt, and run notes carrying the array (or single) job id.

Not this skill: writing or repairing a plugin (the resource specs under
`hvantk/skills/<provider>/SKILL.md`); triaging a drift signal or a broken probe
(`drift-triage`).

## 2. Preflight

1. Site profile:
   `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`;
   stop and ask for any missing value (`skills/workflow-conventions/SKILL.md` section 4).
   Every later command that reads a profile value starts with the same
   `source "$HVANTK_SITE_PROFILE" &&`.
2. Checkout and its commit: `git -C "$HVANTK_REPO" rev-parse HEAD` and
   `git -C "$HVANTK_REPO" status --porcelain`; compare with the commit the user expects.
3. Image and its digest: `readlink -f "$HVANTK_SIF"` and
   `sha256sum "$(readlink -f "$HVANTK_SIF")"`.
4. Inputs and outputs: writable from compute nodes, under `$HVANTK_WORK` or similar, never
   home.
5. What each dataset needs: there is no single command for this.
   `hvantk plugins describe <provider>` prints the name, version, and each dataset's
   domain/backend, skill path and committed fingerprint path -- nothing about artifact
   type, lifecycle or acquisition. Read the manifest itself,
   `hvantk/skills/<provider>/plugin.yaml`: `artifact_type` decides the output extension
   (`.ht` table, `.mt` variant matrix, `.h5ad` expression -- always a single file, never a
   directory -- `.parquet` a pandas-backed table), `lifecycle.parse` (if present) means a
   parse stage and an `--intermediate` path, and `acquisition.mode: byo` means you supply
   the raw files yourself (the reason and any instructions are printed by `reprocess`'s
   own error at build time; see Traps).
6. Egress: whatever the site profile's `HVANTK_EGRESS` says, downloads happen on the login
   node (`references/raw-dir-conventions.md`); the build jobs themselves always run with
   `--skip-download`, so egress changes who runs the download command and when, not the
   procedure below.

## 3. Procedure

1. **Lay out one raw directory per dataset**: `$HVANTK_WORK/raw/<provider>[/<dataset>]/`,
   holding only that dataset's files. Full rules, including dated files and
   `latest_dated`, in `references/raw-dir-conventions.md`.
2. **Stage on the login node, through the image**, uniformly for any dataset with a
   downloader -- `--skip-build` skips only the build; it does not skip the parse stage
   if the dataset declares one, so this also produces the parsed intermediate:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk reprocess provider:dataset --raw-dir "$HVANTK_WORK/raw/provider" \
       --output "$HVANTK_WORK/raw/provider/unused.ht" --skip-build --no-check-drift
   ```
   A download-only `--plugin-arg` (e.g. `genome_build=GRCh38` for `clinvar:variants`'s
   downloader) goes on this call; it is not guaranteed to also be accepted by the
   builder -- a builder's parameters are not a superset of its downloader's (clinvar's
   builder takes `reference_genome`, a different name for the same thing, and has no
   catch-all for an argument it does not declare). For a parse-stage dataset
   (`uniprot-ptm:sites`), give `--intermediate` an explicit path outside the raw
   directory here, so the task line below can reuse it with `--skip-parse
   --intermediate <that same path>` instead of re-parsing inside the job. Then, before
   spending a job on the build, check upstream has not drifted against the committed
   fingerprint: `hvantk drift provider:dataset` through the same `run-in-image.sh`
   line. A `drifted` result whose diff touches `headers` or `checksums` is a schema
   signal: stop and use `drift-triage` first. A clean or `routine` result means it is
   safe to proceed; so does a `probe_failed` caused by a network error ("unavailable,
   proceed"), but any OTHER `probe_failed` (the probe raised on content it did fetch)
   is itself a signal, same as `drifted` -- use `drift-triage` (Traps).
3. **Write the task list**: a TSV, one build per line, absolute paths —
   `<provider:dataset>\t<raw dir>\t<output path>[\t<extra reprocess options>]`. Extra
   options are build-stage ones, split on plain whitespace exactly as written, with no
   quote handling at all (e.g. `--plugin-arg reference_genome=GRCh38`, or
   `--skip-parse --intermediate /abs/path` for a dataset staged with a parse step); a
   value containing a space cannot be expressed in this column.
4. **Submit the array**:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/dataset-build-at-scale/templates/submit-build-array.sh" \
     --list "$HVANTK_WORK/out/tasks.tsv" --name <run>-build --cpus 4 --mem-gb 16 \
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
     --name <run>-verify --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/dataset-build-at-scale/templates/verify-builds.py" \
       --list "$HVANTK_WORK/out/tasks.tsv" --receipt "$HVANTK_WORK/out/verify.json"
   ```
   Details in `references/verification.md`.
7. **Finish the run notes.** If a build refreshed a dataset whose rebuild-ledger row was
   pending (`run-in-image.sh -- hvantk drift --ledger`, same as every other login-node
   command here), run it the same way:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk drift --mark-rebuilt provider:dataset
   ```
   and commit the ledger alongside the change that motivated the rebuild.

## 4. Scale-up

`references/sizing-by-dataset.md` has the full table. In short: one request sizes every
task in the array, so size it for the list's largest dataset — a dataset an
order of magnitude bigger than the rest (dbNSFP next to HGNC, say) belongs in its own
`hpc-hail-job` submission, not in the array. `--max-concurrent` caps how much of the
node share the array can claim at once, independent of how many lines are in the list.
Each task gets its own node-local scratch and its own Hail log directory from
`hail-job.sbatch`, so concurrent tasks never collide even on the same node. Outputs on
the parallel filesystem are directories of many part files; keep each dataset's dated
raw file and its `.provenance.json` sidecar discoverable next to its output.
Re-running only the failed indices (`--tasks`) is the normal way to recover from a
partial array failure, under a fresh job id rather than a retry of the old one.

## 5. Traps

- `expected exactly one file matching ...`. Cause: two files in one raw dir (a stray
  second download). Fix: one dataset per raw directory; give `--intermediate` a path
  outside the raw dir for a parse-stage dataset.
- An HTTP or connection error inside a build job. Cause: compute nodes have no egress.
  Fix: stage the raw file on the login node first; build jobs run with
  `--skip-download`.
- `... requires data you supply yourself (acquisition.mode=byo, reason=...)`. Cause: a
  BYO dataset with an empty or missing `--raw-dir`. Fix: place the raw file yourself;
  this error also names the reason and any instructions the manifest carries.
- `... has no lifecycle.download declared`. Cause: a *different* dataset state from
  BYO -- no `acquisition` block at all, and no downloader written yet (automatable in
  principle, just not done). Fix: `--skip-download` and place the file yourself, or
  write the downloader.
- The output-extension error (`<dataset> writes .ht ...; got --output ...`). Cause: the
  artifact type decides the extension (`.ht` table, `.mt` variant matrix, `.h5ad`
  expression, `.parquet` a pandas-backed table), and a mismatch is rejected before the
  build runs. Fix: match `--output`'s extension to the manifest's `artifact_type`.
- A `--plugin-arg` value arriving as the wrong type. Cause: the coercion rules
  (`true`/`false` to bool, digits to int, a comma to a list) — a single value that
  itself contains a comma cannot be expressed this way. Fix: read the builder's or
  downloader's own parameter docstring for the expected type.
- A build job's log carries a WARNING about the drift probe, and the output's sidecar
  has `source_fingerprint: "probe-unavailable:source-unreachable"`. Cause: the
  platform runs the drift probe on every build regardless of `--no-check-drift` (that
  flag only gates `reprocess`'s own post-build `drift: <status>` print), and a compute
  node has no egress to reach upstream with. Fix: nothing -- the build is fine. The
  real signal is step 2's login-node `hvantk drift provider:dataset`, run *before* the
  build, which compares what upstream serves right now against the committed
  fingerprint (never your staged file); a `probe_failed` there from a network error
  just means "unavailable, proceed" -- but any OTHER `probe_failed` (the probe raised
  on content it did fetch) is itself a signal, same as a `drifted` schema result: use
  `drift-triage`.
- `sacct` showing one row `<A>_[2-3]` for a while after submission. Cause: tasks still
  pending can be reported collapsed into one range row. Fix: wait for every task to
  reach a terminal state, not for that row to disappear.
- One failed task does not fail the array. Cause: Slurm array tasks are independent
  jobs sharing a job id. Fix: check per-task states (`sacct -j <A> -X -n -o
  JobID,State`), never the submission's own exit code alone.
- A built table's values differing from the committed `sample_rows.json`, or an extra
  column the committed `schema.json` does not list. Cause: the committed snapshot is
  fixture-derived, a floor, not an equality target (GenCC's full build carries 17 more
  columns than its 14-column fixture). Fix: `verify-builds.py` checks the snapshot as a
  contract (every declared field present, extra allowed) and the sample keys only
  informationally; compare row values across your own builds, never against the
  fixture's.
- A relative path in the task list. Cause: an array task's working directory is the log
  directory, not the checkout. Fix: absolute paths only, both columns.
- `hvantk: command not found` outside the image, or a command that silently wrote
  into the wrong checkout. Cause and fix: see `workflow-conventions`'s Runtime rules
  (section 4) -- route it through `run-in-image.sh`, never a bare `hvantk` or
  `python -m hvantk`.
- An output moved without its sidecar. Cause: `.provenance.json` sits next to the
  output, not inside it, and does not record the raw file's path (see
  `references/verification.md` for its actual fields). Fix: move both together, or
  neither.
- A task list with a stale raw directory. Cause: yesterday's dated file next to today's
  is tolerated for `latest_dated` plugins (ClinGen, GenCC) but not for one without it
  (HGNC). Fix: keep exactly one matching file per raw directory
  (`references/raw-dir-conventions.md`).

## 6. Verification

1. `sacct -j <A> -X -n -o JobID,State` — every task row `COMPLETED`; `<A>` and `<A>_<i>`
   both count as the array's own ids.
2. `verify-builds.py`'s receipt (`references/verification.md`): per dataset,
   `schema_contract_ok` true, `schema_id_matches` true, the row count and the
   committed snapshot keys' presence (`snapshot_keys_present`; no sample rows)
   recorded for your own comparison against the next build, and `provenance_sha256`
   tying the verdict to exactly this output.
3. Each output carries `_SUCCESS` (or, for a single-file `.h5ad`/`.parquet` artifact,
   exists as a file) and a `.provenance.json` sidecar naming the `schema_id` you
   expected.
4. Run notes carry `commit:`, `container:`, and one `job:` line per job submitted (the
   array counts as one `job:` line; the verification job is another).
