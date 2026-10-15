# Verifying a build

`templates/verify-builds.py` is the one command that answers "did this build actually
produce what the plugin promises", for one dataset or for a whole task list. Run it
inside the image (it needs Hail), as a job:

```bash
source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
  --name <run>-verify --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
  python "$HVANTK_REPO/skills/dataset-build-at-scale/templates/verify-builds.py" \
    --list "$HVANTK_WORK/out/tasks.tsv" --receipt "$HVANTK_WORK/out/verify.json"
```

## What the receipt holds

One line on stdout per dataset (`OK`/`FAIL`, row count, `schema_id`, whether the schema
matched, whether every committed sample key was present), and a receipt JSON with one
entry per dataset under `"datasets"`:

| field | meaning |
| --- | --- |
| `output` | the path verified |
| `success_marker` | whether `_SUCCESS` (or, for a single-file artifact, the file itself) exists |
| `sidecar` | whether `<output>.provenance.json` exists |
| `schema_id` | the sidecar's own `schema_id` |
| `provenance_sha256` | SHA-256 of the sidecar's bytes |
| `schema_matches` | the output's live schema equals the plugin's committed `schema_snapshot` |
| `row_count` | rows (Hail) or `n_obs` (AnnData) |
| `sample_keys_present` | every key in the plugin's committed `row_snapshot` resolved to a row |
| `sample_rows` | the collected rows for those keys (values, not just presence) |

The command exits 1 if any dataset is missing, lacks `_SUCCESS`, fails its schema
comparison, or is missing a committed key. A verdict about one run is bound to its
`provenance_sha256`: if you rebuild the same dataset later and the sidecar's hash
changes, this receipt is about the earlier build, not the one on disk now.

## How the snapshots relate to a full build

A plugin's `tests/snapshots/schema.json` and `tests/snapshots/sample_rows.json` are
generated from its tiny committed test fixture (`hvantk/skills/<provider>/tests/testdata/`),
not from a full upstream file. The **schema** is stable between the fixture and a full
build — the builder's column renaming and typing do not depend on row count — so
`schema_matches` is a real, exact check. The **sample-row values** the snapshot carries
can differ for a full build of the same keys: dates roll forward, an array field that
was short in the fixture may be wider in the real file. That is why `verify-builds.py`
only checks that the committed sample *keys* are present (`sample_keys_present`), and
reports the collected row *values* in the receipt for you to look at, rather than
asserting they equal the fixture snapshot's values. A key that is missing after a
rebuild is still a real problem: it means the upstream file no longer carries a row the
plugin's own tests depend on, or the builder's renaming silently dropped it (HGNC logs a
missing rename at WARNING and keeps going unless `strict_fields=true` — exit 0 is not
proof the schema didn't change).

## Comparing two builds of the same dataset

To check that a rebuild reproduces an earlier one (after an upstream refresh you expect
to be a no-op, or after a plugin change you expect not to affect this dataset), compare
three things, in order of how much they tell you:

1. **Row count.** The cheapest check; a changed count means upstream added or removed
   rows, or a filter in the builder changed.
2. **Sample-row values**, for the same keys, byte for byte. `verify-builds.py`'s
   `sample_rows` in the receipt is already in a JSON-stable form
   (`hvantk.tests._snapshot_utils.collect_sample_rows`'s output); diff two receipts'
   `sample_rows` for the same dataset directly.
3. **Schema.** The least likely to change day to day; a difference here is the one worth
   treating as news rather than noise (a renamed or retyped column upstream).

## The provenance sidecar is what a verdict binds to

`<output>.provenance.json` carries `plugin`, `dataset`, `plugin_version`,
`source_fingerprint`, `schema_id`, `build_timestamp` and `builder_commit`. It is the
only place that records which raw file and which commit produced a given output
directory — not the output directory's own contents, and not its mtime. Keep the
sidecar next to its output when you move, copy or archive a build; a verification
receipt or a drift-triage rebuild record that cites a `provenance_sha256` is citing this
file, not the table.
