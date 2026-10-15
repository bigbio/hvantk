# Classifying a drift signal

The decision table this skill's step 2 uses. See `hvantk/skills/_conventions/SKILL.md`
section 12 for the probe contract and the automated drift workflow this table is read
against; this file is the operational read of it, not a restatement.

## Reading a report

`hvantk drift <provider:dataset> --json` (or the offline form's `drift ... -- --json`)
prints a list with one entry for the dataset. The fields that matter:

| field | meaning |
| --- | --- |
| `status` | `clean`, `drifted`, `probe_failed`, or `stub` |
| `diff.added` / `diff.removed` | top-level fingerprint keys gained or lost |
| `diff.changed` | keys present on both sides with a different value |
| `probe_error` | set only for `probe_failed`; the exception the probe raised |

A drift-bot pull request's body carries the same `diff` rendered as a markdown table
before the raw JSON; it is derived from this report, not a separate signal.

## Signal -> action -> bookkeeping

| signal | what it means | action | bookkeeping |
| --- | --- | --- | --- |
| `diff.changed` holds `headers` or `checksums`, or `diff.added`/`diff.removed` is non-empty | a **schema** change: the column list or header hash moved, or a top-level key appeared or disappeared | read the plugin's rename map and `references/repair-patterns.md`; repair before accepting | `--regenerate --accepted-in "PR #N"` once the repair is validated; ledger row `signal: schema` |
| only `extras` / `source_version` / `informational` moved (everything else equal once `fetched_at`/`probe_version`/`informational` are stripped) | **routine**: content changed, schema did not | usually nothing to repair; rebuild if the dataset is due for a refresh | `--regenerate --accepted-in "PR #N"`; ledger row `signal: routine` |
| `status: probe_failed`, `probe_error` set | the probe itself is broken (network, parsing, a missing baseline) -- **not** a drift signal | fix `drift_probe.py`; bump `PROBE_VERSION` if its output shape changed | none for the bump alone -- `record_regeneration` writes no row for a `probe_version`-only change |
| `status: stub` | a documentation-only source with no programmatic probe | nothing to do | none |

Three outcomes, not two: accept as-is (routine, or a schema change nothing downstream
reads yet), repair then accept (a schema change that broke something), or fix the probe
(the signal was never about the dataset at all).

## What the drift bot does and does not do

`.github/workflows/drift.yml` runs `hvantk drift --all --json` on a schedule and opens
one pull request per signal: every routine signal in a run goes onto one shared
`drift/routine-batch` branch and PR; a schema signal keeps its own branch and PR. Across
every drift PR the bot has opened, it has only ever bumped a fingerprint and recorded a
ledger row -- it has never touched a builder, a manifest, a probe, or a snapshot.
Merging a schema-signal PR green therefore accepts the new fingerprint; it does **not**
mean the plugin still builds. Run this skill's repair and rebuild steps regardless of
which PR merged the fingerprint, whenever the signal is a schema one.

`probe_failed` and `stub` entries are logged by the bot but never turned into a pull
request -- the former is an infrastructure failure, the latter has no programmatic
signal to act on. Seeing either reported by `hvantk drift` directly, outside the bot,
means the same thing: it is not a dataset-content question.
