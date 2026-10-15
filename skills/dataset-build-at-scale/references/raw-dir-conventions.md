# Raw-directory conventions

`hvantk reprocess` hands a no-parse dataset's builder the raw directory itself, and the
builder resolves it to one file with `hvantk/core/utils/raw_inputs.py::resolve_raw_file`.
Get the directory wrong and the error surfaces inside the build, not before it.

## One file per dataset directory

A raw directory holds exactly one dataset's raw file(s). Never share a directory
between two datasets, even two from the same provider: `resolve_raw_file` matches
`fnmatch` patterns against every non-hidden regular file in the directory (no
recursion), and a second matching file turns "found 1" into "found 2" — the same error
whether the second file is a stray download or yesterday's run.

Lay raw directories out as `$HVANTK_WORK/raw/<provider>[/<dataset>]/`: one subdirectory
per dataset for a multi-dataset provider (its datasets' patterns can otherwise collide),
a bare `$HVANTK_WORK/raw/<provider>/` for a single-dataset one.

## Patterns are per plugin

Each plugin's builder calls `resolve_raw_file` with the patterns that match what its own
downloader writes (and what its test fixture is named), so do not assume a provider's
raw directory accepts just anything: `*.txt`/`*.tsv` for HGNC, `*.csv` for ClinGen,
`*.tsv` for GenCC, each case-insensitive. Hidden files (`.`-prefixed) are always
ignored, so a stray `.DS_Store` never counts as a second candidate. A *visible*
leftover can still appear, though: `hvantk/core/utils/http.py::stream_to_file` writes
to a temporary `<name>.<random>` path and renames it into place only on success,
cleaning the temp file up itself on any raised exception (an HTTP error, a short
read) -- but an ungraceful interruption (the process killed outright, not an
exception) skips that cleanup and leaves the temp file sitting in the raw dir as a
second, non-hidden match.

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

A dataset whose manifest declares `lifecycle.parse` (`uniprot-ptm:sites` is one) reads
its raw directory at the parse stage and writes a parsed intermediate. `--intermediate`
is optional: omitted, it defaults to `<raw-dir>/<provider>_<dataset>.intermediate`,
which is *inside* the raw directory — exactly the kind of second file
`resolve_raw_file` cannot tell from a stray download on the next run. Give
`--intermediate` an explicit path outside the raw directory for any dataset you will
`reprocess` more than once from the same raw directory, which in practice means every
parse-stage dataset in this skill's procedure: stage it once with that explicit path,
then reuse the same path from the task list with `--skip-parse --intermediate <path>`
so the array task consumes it directly instead of re-parsing.

## BYO datasets and `acquisition`

A dataset with `acquisition.mode: byo` has no `lifecycle.download`; `reprocess` skips
the download stage implicitly (passing `--skip-download` is also legal) but still
requires `--raw-dir` to be a non-empty directory unless you also pass `--skip-parse`
with an explicit `--intermediate` (then the raw dir is never opened). You supply the raw
file yourself, named so the plugin's patterns still match it. There is no single
command that prints a BYO dataset's reason or instructions ahead of time: they appear
in `reprocess`'s own error when `--raw-dir` is empty (see the skill's Traps) -- read
`hvantk/skills/<provider>/plugin.yaml`'s `acquisition:` block directly to find out
before that point.

## `--plugin-arg` forwarding is not stage-agnostic

Every `--plugin-arg KEY=VALUE` on a `reprocess` invocation is forwarded, as a coerced
kwarg, to whichever of download, parse and build actually run -- but a stage's function
only accepts the parameters it declares, and a builder commonly has no catch-all for an
argument it does not recognize. `clinvar`'s download function takes `genome_build`
(and has `**kwargs`, so an unrecognized extra is harmless there); its builder takes
`reference_genome` -- a different name for the same thing -- with no `**kwargs` at all,
so a `--plugin-arg genome_build=...` aimed at the download stage would raise a
`TypeError` if it reached the builder. Do not assume one `--plugin-arg` list is safe
everywhere: pass staging-only arguments on the `--skip-build` staging call (step 2) and
build-only arguments in the task list's own extra-options column (step 3), not both in
the same place.

## Where downloads happen

Compute nodes on a typical site have no internet access (`HVANTK_EGRESS` in the site
profile says so). Every download therefore happens on the login node, through the
image, before any job is submitted:

```bash
source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
  hvantk download hgnc --output "$HVANTK_WORK/raw/hgnc/hgnc_complete_set.txt"
```

The array job itself always runs with `--skip-download --no-check-drift`: by the time a
build job starts, the raw file (and, for a parse-stage dataset, the intermediate) must
already be staged.

## Refreshing a stale raw file

A dataset's downloader refuses to overwrite an existing file; add
`--plugin-arg overwrite=true` to replace it with a fresh copy in place (or, for
`clingen`/`gencc`, pass `--overwrite` to the plugin's own `hvantk download <provider>`
command) -- no manual removal needed first, for any of the three. Refreshing leaves the
old dated file in place for `latest_dated` plugins, which is fine (see above); HGNC has
only the one filename, so `overwrite=true` simply replaces it.

## The drift check before a build

`hvantk drift <provider:dataset>` on the login node (through `run-in-image.sh`)
compares what upstream serves *right now* against the dataset's committed fingerprint
-- it reaches out to the real source itself; it never reads your staged file, so it can
tell you upstream has moved even before you have re-staged anything. A `schema` signal
(a column added, removed or renamed in the file's header) means the plugin itself may
need a repair before the build is worth trusting — see `drift-triage`. A `routine`
signal, or a clean result, means it is safe to proceed. `probe_failed` (no egress from
the login node either, or a transient error) means the signal is unavailable, not that
anything is wrong: proceed, and run the check again later when connectivity allows it.
