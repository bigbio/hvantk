# Verifying a build

`templates/verify-builds.py` is the one command that answers "did this build actually
satisfy what the plugin promises", for one dataset or for a whole task list. Run it
inside the image (it needs Hail), as a job:

```bash
source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
  --name <run>-verify --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
  python "$HVANTK_REPO/skills/dataset-build-at-scale/templates/verify-builds.py" \
    --list "$HVANTK_WORK/out/tasks.tsv" --receipt "$HVANTK_WORK/out/verify.json"
```

## What "satisfy" means: a contract, not equality

A plugin's committed `tests/snapshots/schema.json` is a floor, not a target to match
exactly: it is generated from a tiny committed test fixture
(`hvantk/skills/<provider>/tests/testdata/`), and a full build from the real upstream
file routinely carries more columns than the fixture did. GenCC's upstream file has 31
columns; its committed snapshot declares 14 row fields; a full build legitimately has
17 more fields than the snapshot lists. `verify-builds.py` therefore checks that every
field the snapshot declares is present in the built schema with the same type
(`schema_contract_ok`) and separately lists whatever extra fields the build has
(`extra_fields`) -- an extra field is not a failure, a missing or retyped one is.

The committed `tests/snapshots/sample_rows.json`'s *keys* are informational only
(`snapshot_keys_present`), never a pass/fail condition: some plugins' fixtures are
hand-authored rather than sub-sampled from a real file, so their keys can be
identifiers that do not exist anywhere upstream. ClinGen's fixture invents its own
MONDO ids for exactly this reason: a build from the real upstream file is not
guaranteed to contain those same gene-disease pairings, so treating a missing key as a
hard failure would punish the mismatch between a hand-authored fixture and reality, not
an actual problem with the build. A missing key here tells you the fixture and the
real file disagree on identifiers, not that the build is wrong.

## What the receipt holds

One line on stdout per dataset (`OK`/`FAIL`, row count, `schema_id`, whether it matches
the registry's, whether the schema contract holds, the extra fields, whether the
snapshot's sample keys were found), and a receipt JSON with one entry per dataset under
`"datasets"`, keyed `"<dataset>=<output>"` (so verifying the same dataset's output
twice, e.g. before and after a rebuild, does not collide):

| field | meaning |
| --- | --- |
| `output` | the path verified |
| `success_marker` | whether `_SUCCESS` (or, for a single-file `.h5ad`/`.parquet` artifact, the file itself) exists |
| `sidecar` | whether `<output>.provenance.json` exists |
| `schema_id` | the sidecar's own `schema_id` |
| `schema_id_matches` | whether it equals the registry's `schema_id` for this dataset |
| `provenance_sha256` | SHA-256 of the sidecar's bytes |
| `schema_contract_ok` | every field the committed snapshot declares is present with the same type |
| `missing_fields`, `type_mismatches`, `extra_fields` | the contract check's detail |
| `row_count` | rows (Hail), `n_obs` (AnnData) or row count (a parquet table) |
| `snapshot_keys_present` | `{"found": n, "total": m}` for the committed `row_snapshot`'s keys -- informational, see above |

One dataset's exception never costs the others' results: each is wrapped in its own
try/except, so a dataset that errors still leaves a receipt entry (with an `error`
field) and the others are still verified. The command exits 1 if any dataset's output
is missing, lacks `_SUCCESS`/its sidecar, fails the schema contract, or has a
`schema_id` that does not match the registry's -- never on the informational key
check alone. A verdict about one run is bound to its `provenance_sha256`: if you
rebuild the same dataset later and the sidecar's hash changes, this receipt is about
the earlier build, not the one on disk now.

## Comparing two builds of the same dataset

To check that a rebuild reproduces an earlier one (after an upstream refresh you expect
to be a no-op, or after a plugin change you expect not to affect this dataset), compare
three things, in order of how much they tell you, against your OWN earlier receipt or
build -- never against the fixture's committed `sample_rows.json`, whose values (and,
for a hand-authored fixture, whose very identifiers) are not guaranteed to resemble a
full build's:

1. **Row count.** The cheapest check; a changed count means upstream added or removed
   rows, or a filter in the builder changed.
2. **The schema contract's own fields**: `extra_fields` changing between two of your
   own builds (not relative to the snapshot, which is expected to differ) is worth
   investigating; `missing_fields`/`type_mismatches` appearing at all is a regression.
3. **Sample-row values**, for keys you choose yourself from the live table (the
   snapshot's own sample keys may not even resolve, per ClinGen above).

## The provenance sidecar is what a verdict binds to

`<output>.provenance.json` carries `plugin`, `dataset`, `plugin_version`,
`source_fingerprint`, `schema_id`, `build_timestamp`, `builder_commit` (the hvantk code
that built it: the commit of the checkout a task ran, which `hail-job.sbatch` passes into
the image as `HVANTK_BUILDER_COMMIT`, with `-dirty` for uncommitted edits under
`hvantk/`; for a build by the image's own copy of hvantk, the commit the image was built
from; else `hvantk <version>`, as with a wheel install or a task on a node where git
cannot answer), `parents` and `build_parameters`. It does **not**
record the raw file's path: `source_fingerprint` is the drift probe's fingerprint (or,
built through a job with no egress, the `probe-unavailable:source-unreachable` token --
see the skill's Traps), not a path, and nothing else on the sidecar names one either.
Keep the sidecar next to its output when you move, copy or archive a build; a
verification receipt or a drift-triage rebuild record that cites a `provenance_sha256`
is citing this file, not the table.
