# Recording a drift decision

The fingerprint and the rebuild ledger, and when the offline form stands in for the
real probe. `hvantk/skills/_conventions/SKILL.md` section 12 is the authoritative
contract (the probe's return shape, the keys a comparison ignores, the two risk tiers,
the automated workflow); this file is the operational sequence, not a restatement.

## The two files

- **Fingerprint**: `hvantk/skills/<provider>/[<dataset>/]tests/drift_fingerprint.json`
  (clingen and hgnc each ship one at `tests/drift_fingerprint.json`, one dataset per
  provider). The committed baseline `hvantk drift` compares the live probe's output
  against.
- **Ledger**: `hvantk/resources/drift_ledger.json`, one row per dataset:
  `last_upstream_change`, `accepted_in`, `signal` (`routine` or `schema`), `rebuilt_at`
  (null until marked rebuilt).

Both files are written only through `hvantk.core.plugin.drift_runner` /
`drift_ledger`, never by hand: a hand-edited fingerprint is rejected as a baseline that
was never captured from a live probe (and so can never equal a real observation again),
and a hand-edited ledger row is just as easy to get silently wrong as any other
hand-maintained provenance record.

## Sequence

1. `hvantk drift --regenerate <provider:dataset> --accepted-in "<PR or branch>"` --
   overwrites the fingerprint with the live probe's current output, and records a
   ledger row IF the new baseline differs materially from the old one (ignoring
   `fetched_at`/`probe_version`/`informational`). No row is written for a first capture
   or for a `probe_version`-only change; add one by hand in that case if upstream
   genuinely also moved in the same window. `--accepted-in` defaults to `manual` --
   always pass the real PR or branch instead. `--no-ledger-row` skips the row (the
   drift bot regenerates this way and records its own row at its own point once its
   checks pass; a human triaging by hand should not normally need it).
2. `hvantk drift --mark-rebuilt <provider:dataset>` -- sets the row's `rebuilt_at` to
   now. Refuses a dataset with no ledger row: there is nothing pending to clear.
3. `hvantk drift --ledger` -- lists every dataset never rebuilt, or rebuilt before its
   own `last_upstream_change`. The dataset just handled must not appear in that list.

Run 1 and 2 in that order: `--mark-rebuilt` requires the row `--regenerate` writes, and
fails outright if it runs first.

## The offline form

`templates/probe-offline.py drift --dataset <provider:dataset> --file <path> --
--regenerate --accepted-in "<PR or branch>"` does the same as step 1, except the probe's
`requests.head`/`requests.get` calls are answered from `<path>` instead of the network.
Reach for it when:

- the machine running the triage has no egress (a compute node, a sandboxed
  evaluation), or
- the baseline should be pinned to the EXACT file the rebuilt dataset was validated
  against, rather than whatever upstream happens to serve at the moment `--regenerate`
  is run -- upstream can change again in the gap between the drift report and the
  regenerate call, and the committed baseline should describe the file the dataset
  actually came from, not a third, unvalidated fetch.

Coverage is probes that call the module-level `requests.head`/`requests.get` functions
directly -- clingen and hgnc both do. A probe built on `requests.Session` or any other
HTTP client is not covered; `templates/probe-offline.py`'s own module docstring says so.

## What the pull request must carry

The diff should contain, and no more than: the repair (if any), the regenerated
`drift_fingerprint.json`, the ledger row, and the fixture or snapshots ONLY if they were
regenerated on purpose -- stated as such in the commit message, per
`references/repair-patterns.md`. Subject line `fix(<provider>): accept <what changed
upstream>`; the body says what moved upstream and, if anything was repaired, what the
repair was.
