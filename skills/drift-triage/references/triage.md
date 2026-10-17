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

`observed` and `expected` are printed in full alongside `diff`, including
`fetched_at`/`probe_version`/`informational` -- but those three keys are stripped out
BEFORE `diff` is computed, so `informational` (where a probe records something like
`last_modified` for a human reader) can never appear in `diff` at all, under any
status. Comparing the raw `observed`/`expected` blocks by eye and reacting to a
difference there, rather than reading `diff.changed`, is the trap this causes.

A drift-bot pull request's body carries the same `diff` rendered as a markdown table
before the raw JSON; it is derived from this report, not a separate signal.

## Signal -> action -> bookkeeping

| signal | what it means | action | bookkeeping |
| --- | --- | --- | --- |
| `diff.changed` holds `headers` or `checksums`, or `diff.added`/`diff.removed` is non-empty | a **schema** change: the column list or header hash moved, or a top-level key appeared or disappeared | read the plugin's rename map and `references/repair-patterns.md`; repair before accepting | `--regenerate --accepted-in "PR #N"` once the repair is validated; ledger row `signal: schema` |
| only `extras` / `source_version` moved (everything else equal once `fetched_at`/`probe_version`/`informational` are stripped) | **routine**: content changed, schema did not | usually nothing to repair; rebuild if the dataset is due for a refresh | `--regenerate --accepted-in "PR #N"`; ledger row `signal: routine` |
| `status: probe_failed`, `probe_error` set | the probe itself is broken -- **not** a drift signal; see "probe_failed, by cause" below | depends on the cause | none for a cause that is really a probe change -- `record_regeneration` writes no row for a `probe_version`-only change |
| `status: stub` | a documentation-only source with no programmatic probe | nothing to do | none |

Three outcomes for a schema signal, not two: accept as-is (nothing downstream reads the
changed column yet), repair then accept (something broke), or -- if the signal turns
out not to be about the dataset at all -- fix the probe instead.

## `probe_failed`, by cause

`probe_failed` is one status covering three different problems, and the fix is
different for each. `probe_error` is the only field that distinguishes them, so read
its text before doing anything:

| cause | what `probe_error` looks like | fix |
| --- | --- | --- |
| no egress from where the probe ran | an HTTP / connection error | run from a login node with egress, or answer `templates/probe-offline.py` from a file already downloaded |
| no committed baseline, or one that was never captured from a live probe | "missing expected fingerprint" or "... was never captured from a live probe" | `--regenerate` |
| the probe's own code is broken | a parse error, a missing marker, anything that is not network- or baseline-shaped | fix `hvantk/skills/<provider>/drift_probe.py`; bump `PROBE_VERSION` only if the probe's OUTPUT shape changed -- any version difference at all suppresses the next `--regenerate`'s ledger row (`references/bookkeeping.md`) |

A network error or a missing/never-captured baseline (the first two causes) is not
itself a signal about the dataset -- proceed once egress or a baseline exists. The
third cause is different: a probe that raised while parsing content it DID
successfully fetch (a missing header marker, say) is not something to proceed past
the way those two are. It IS a signal, same standing as a `drifted` schema result --
triage it, don't wave it off as "just an infra hiccup".

## What the drift bot does and does not do

`.github/workflows/drift.yml` runs `hvantk drift --all --json` on a schedule and opens
one pull request per signal: every routine signal in a run goes onto one shared
`drift/routine-batch` branch and PR; a schema signal keeps its own branch and PR. The
bot only ever bumps a fingerprint and records a ledger row -- it never touches a
builder, a manifest, a probe, or a snapshot. Merging a schema-signal PR green
therefore accepts the new fingerprint; it does **not** mean the plugin still builds.
Run this skill's repair and rebuild steps regardless of which PR merged the
fingerprint, whenever the signal is a schema one.

`probe_failed` and `stub` entries are logged by the bot but never turned into a pull
request -- the former is an infrastructure failure, the latter has no programmatic
signal to act on.
