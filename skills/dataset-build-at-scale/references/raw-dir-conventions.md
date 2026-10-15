# Raw-directory conventions

`hvantk reprocess` hands a no-parse dataset's builder the raw directory itself, and the
builder resolves it to one file with `hvantk/core/utils/raw_inputs.py::resolve_raw_file`.
Get the directory wrong and the error surfaces inside the build, not before it.

## One file per dataset directory

A raw directory holds exactly one dataset's raw file(s). Never share a directory
between two datasets, even two from the same provider: `resolve_raw_file` matches
`fnmatch` patterns against every non-hidden regular file in the directory (no
recursion), and a second matching file turns "found 1" into "found 2" — the same error
whether the second file is a stray download, a `.intermediate` file from a parse stage,
or yesterday's run.

Lay raw directories out as `$HVANTK_WORK/raw/<provider>[/<dataset>]/`: one subdirectory
per dataset for a multi-dataset provider (its datasets' patterns can otherwise collide),
a bare `$HVANTK_WORK/raw/<provider>/` for a single-dataset one.

## Patterns are per plugin

Each plugin's builder calls `resolve_raw_file` with the patterns that match what its own
downloader writes (and what its test fixture is named), so do not assume a provider's
raw directory accepts just anything: `*.txt`/`*.tsv` for HGNC, `*.csv` for ClinGen,
`*.tsv` for GenCC, each case-insensitive. Hidden files (`.`-prefixed) are always
ignored, so a stray `.DS_Store` or a partial `.tmp` download never counts as a second
candidate.

## Dated files and `latest_dated`

ClinGen and GenCC publish real-time snapshots: their downloaders write
`<prefix>-<YYYY-MM-DD>.<ext>`, and their plugins set `latest_dated=True`. With more than
one file matching the patterns, `resolve_raw_file` only resolves to one of them when
every match is the same file under a different date (the stem with the date string
removed is identical across all of them) and the latest date is not itself a tie; it
then uses the newest and logs the others as ignored at WARNING. A plugin without
`latest_dated` (HGNC) still raises on a second file, dated or not — don't leave
yesterday's HGNC file next to today's.

An unrelated dated file that happens to match the same patterns (a different report
also named `<other>-<date>.tsv`) breaks the "same file, different date" check and falls
through to the same "found N" error — which is correct: there is no uniform rule for
picking between two different things.

## Parse-stage plugins and `--intermediate`

A dataset whose manifest declares `lifecycle.parse` reads its raw directory at the parse
stage and writes a parsed intermediate to `--intermediate` (defaulting to
`<raw-dir>/<provider>_<dataset>.intermediate` when omitted, which is *inside* the raw
directory — exactly the kind of second file `resolve_raw_file` cannot tell from a stray
download on the next run). Give `--intermediate` a path outside the raw directory for
any dataset you will `reprocess` more than once from the same raw directory.

## BYO datasets and `acquisition`

A dataset with `acquisition.mode: byo` has no `lifecycle.download`; `reprocess` skips
the download stage implicitly (passing `--skip-download` is also legal) but still
requires `--raw-dir` to be a non-empty directory unless you also pass `--skip-parse`
with an explicit `--intermediate` (then the raw dir is never opened). You supply the raw
file yourself, named so the plugin's patterns still match it; `hvantk plugins describe
<provider>` prints the acquisition reason and any instructions the manifest carries.

## `--plugin-arg` forwarding

Every `--plugin-arg KEY=VALUE` on a `reprocess` invocation is forwarded, as a coerced
kwarg, to whichever of download, parse and build actually run. A builder-only argument
(e.g. `reference_genome=GRCh38`) is harmless to also pass at a download-only stage (the
lifecycle functions accept `**kwargs` and ignore what they don't use) — there is no need
to split one `reprocess` call's extras by stage.

## Where downloads happen

Compute nodes on a typical site have no internet access (`HVANTK_EGRESS` in the site
profile says so). Every download therefore happens on the login node, through the image,
before any job is submitted:

```bash
source "$HVANTK_SITE_PROFILE" && "${HVANTK_CONTAINER_BIN:-singularity}" exec \
  -B "$HVANTK_WORK:$HVANTK_WORK" -B "$HVANTK_REPO:$HVANTK_REPO" --env PYTHONPATH="$HVANTK_REPO" \
  "$HVANTK_SIF" hvantk download hgnc --output "$HVANTK_WORK/raw/hgnc/hgnc_complete_set.txt"
```

The array job itself always runs with `--skip-download --no-check-drift`: by the time a
build job starts, the raw file must already be staged.

## Refreshing a stale raw file

A dataset's downloader refuses to overwrite an existing file; add
`--plugin-arg overwrite=true` to replace it with a fresh copy (or, for `clingen`/`gencc`,
pass `--overwrite` to the plugin's own `hvantk download <provider>` command). Refreshing
leaves the old dated file in place for `latest_dated` plugins — which is fine, see
above — but for HGNC, remove the old file first.

## The drift check before a build

`hvantk drift <provider:dataset>` on the login node (through the image) compares the
freshly staged file against the dataset's committed fingerprint before you spend a job
building it. A `schema` signal (a column added, removed or renamed in the file's header)
means the plugin itself may need a repair before the build is worth trusting — see
`drift-triage`. A `routine` signal, or a clean result, means it is safe to proceed.
