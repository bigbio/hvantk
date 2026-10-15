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

The fingerprint is written only through the tooling, never by hand. What "rejected"
actually means is narrower than it sounds: only a baseline shaped like a placeholder --
carrying the specific markers a hand-seeded file leaves (an empty checksum, a sentinel
string, `fetched_at` at the Unix epoch) -- is refused outright as `probe_failed`. Any
OTHER hand-edit is accepted as a file -- there is no general "this was not written by
the tool" detector -- it just then never equals a live observation again, so the
dataset reports `drifted` on every run forever. Either way the fix is the same:
`--regenerate`, never a manual edit.

The ledger is written the same way in the NORMAL case, but it has one documented,
narrow exception: a hand-added row, for the specific case where `--regenerate` records
none because the two baselines differ only in `probe_version` (step 1 below) while
upstream genuinely also moved in that same window -- the dbnsfp and uniprot-ptm
fingerprint bumps both did this by hand. Outside that one case, write it through
`--regenerate`/the bot, not by hand: a hand-maintained provenance record is exactly as
easy to get silently wrong as any other one.

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

Run 1 BEFORE 2, always -- but the risk is subtler than "2 fails without 1's row". If
the dataset has no ledger row at all yet, running `--mark-rebuilt` first does fail
outright, obviously wrong. The REAL mechanism: `hvantk drift --regenerate` that
records a row (`ledger_update` in `hvantk/core/plugin/drift_ledger.py`) REPLACES the
dataset's row outright -- a fresh `last_upstream_change`, and `rebuilt_at` reset to
`null` -- so a `--mark-rebuilt` run before it is simply overwritten by the regenerate
that follows; wrong order alone is self-correcting, not silently dangerous.

The genuinely false-clean case is different: `--regenerate` records NO row at all (a
first capture, or a `probe_version`-only change -- step 1 above) while
`clingen:gene-disease` (or any dataset that has drifted before) already HAS an older
row from some earlier, unrelated signal. `--mark-rebuilt` then marks THAT stale row
rebuilt, succeeds with no error, and `--ledger` reports clean even though the signal
you are triaging right now was never recorded at all. The only way to notice is to
check that the row's `last_upstream_change` is from THIS run, not a stale one -- which
is also why step 3's `--ledger` passing is necessary but not sufficient; read the row
itself too.

## The offline form

`templates/probe-offline.py drift --dataset <provider:dataset> --file <path> --
--regenerate --accepted-in "<PR or branch>"` does the same as step 1, except the
probe's HTTP calls (see coverage, below) are answered from `<path>` instead of the
network. Reach for it when:

- the machine running the triage has no egress (a compute node, a sandboxed
  evaluation), or
- the baseline should be pinned to the EXACT file the rebuilt dataset was validated
  against, rather than whatever upstream happens to serve at the moment `--regenerate`
  is run -- upstream can change again in the gap between the drift report and the
  regenerate call, and the committed baseline should describe the file the dataset
  actually came from, not a third, unvalidated fetch.

Coverage is probes that call the module-level `requests.head`, `requests.get` or
`requests.request` functions -- clingen and hgnc call `head`/`get` directly; gencc goes
through `hvantk.core.utils.http.request_with_retry`, which calls `requests.request`, so
it is served too. A probe built on its own `requests.Session` (gnomad-metrics) is NOT
served: it is actively BLOCKED for the duration, so it fails loudly as `probe_failed`
rather than silently reaching the real network or recording a fingerprint for the wrong
file. `templates/probe-offline.py`'s own module docstring has the exact mechanism.

## What the pull request must carry

The diff should contain, and no more than: the repair (if any), the regenerated
`drift_fingerprint.json`, the ledger row, and the fixture or snapshots ONLY if they were
regenerated on purpose -- stated as such in the commit message, per
`references/repair-patterns.md`. Subject line `fix(<provider>): accept <what changed
upstream>`; the body says what moved upstream and, if anything was repaired, what the
repair was.
