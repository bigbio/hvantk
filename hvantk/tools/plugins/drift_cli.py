"""`hvantk drift ...` command for fingerprint-based drift detection."""

from __future__ import annotations

import dataclasses
import json
from datetime import datetime
from pathlib import Path

import click

from hvantk.core.plugin import drift_ledger, drift_runner, loader as plugin_loader
from hvantk.core.plugin.api import DriftProbeError, PluginLoadError


EXIT_CLEAN = 0
EXIT_DRIFTED = 1
EXIT_PROBE_FAILED = 2
EXIT_REGISTRY_ERROR = 3
# A caller-side invocation mistake (an invalid flag combination, an unresolvable
# --domain) rather than anything about a dataset's drift status. click.UsageError's
# own exit_code is 2 -- the same as EXIT_PROBE_FAILED -- and both drift.yml and
# drift-health.yml tolerate any `rc <= 2` as "a signal, not a failure" (drifted /
# probe_failed), so a typo'd flag combination used to pass the scheduled workflows'
# own gate instead of failing the run. Raised via _usage_error() below instead of
# click.UsageError so it exits with this code rather than click's.
EXIT_USAGE_ERROR = 4

# hvantk/resources/drift_ledger.json: the rebuild ledger `--regenerate` and the drift bot
# write (see hvantk/core/plugin/drift_ledger.py). Read from this module at call time, so
# a test can point every ledger read and write here at a temporary file.
LEDGER_PATH = drift_ledger.LEDGER_PATH


def _usage_error(message: str) -> None:
    """Report an invalid invocation and exit EXIT_USAGE_ERROR.

    Mirrors what click.UsageError("message") would print -- none of the call sites
    below attach a ctx, so click's own UsageError.show() would skip the "Usage: ..."
    banner anyway and print only this "Error: ..." line -- but exits with its own
    code instead of click's (always 2, colliding with EXIT_PROBE_FAILED; see
    EXIT_USAGE_ERROR above).
    """
    click.echo(f"Error: {message}", err=True)
    raise SystemExit(EXIT_USAGE_ERROR)


@click.command(name="drift")
@click.argument("dataset", required=False)
@click.option(
    "--all", "all_flag", is_flag=True, help="Run drift check for every dataset"
)
@click.option("--domain", default=None, help="Filter by domain (with --all)")
@click.option("--json", "as_json", is_flag=True, help="Emit machine-readable JSON")
@click.option(
    "--regenerate",
    is_flag=True,
    help=(
        "Overwrite drift_fingerprint.json with observed probe output, and record a "
        "rebuild-ledger row if upstream materially changed (not for a first capture "
        "or a probe_version change)"
    ),
)
@click.option(
    "--timeout", default=60, show_default=True, help="Probe timeout in seconds"
)
@click.option(
    "--ledger",
    "ledger_flag",
    is_flag=True,
    help="List datasets whose upstream moved since their last rebuild",
)
@click.option(
    "--mark-rebuilt",
    "mark_rebuilt",
    default=None,
    metavar="DATASET",
    help="Record DATASET as rebuilt now, clearing it from --ledger's stale list",
)
@click.option(
    "--accepted-in",
    "accepted_in",
    default=None,
    metavar="REF",
    help=(
        "With --regenerate: the PR or branch that accepts the new fingerprint, recorded "
        f"in its rebuild-ledger row [default: {drift_ledger.MANUAL_ACCEPTED_IN}]"
    ),
)
@click.option(
    "--no-ledger-row",
    "no_ledger_row",
    is_flag=True,
    help=(
        "With --regenerate: write no rebuild-ledger row (the drift bot records its own "
        "once its checks pass)"
    ),
)
def drift_cmd(
    dataset,
    all_flag,
    domain,
    as_json,
    regenerate,
    timeout,
    ledger_flag,
    mark_rebuilt,
    accepted_in,
    no_ledger_row,
):
    """Compare a plugin's live drift-probe fingerprint against the expected file."""
    # Both only qualify --regenerate's ledger row, so on their own -- or with --ledger /
    # --mark-rebuilt, which never regenerate -- they would be silently ignored.
    if (accepted_in is not None or no_ledger_row) and not regenerate:
        _usage_error("--accepted-in and --no-ledger-row only apply to --regenerate")
    if accepted_in is not None and no_ledger_row:
        _usage_error("--accepted-in and --no-ledger-row are mutually exclusive")
    if accepted_in is not None and not accepted_in.strip():
        _usage_error("--accepted-in needs a PR or branch name")

    if mark_rebuilt is not None:
        # A standalone action, like --regenerate: it takes its OWN dataset name, so a
        # positional `dataset` / --all / --regenerate alongside it would be ambiguous
        # about which dataset is meant. Mutually exclusive with --ledger too -- marking
        # and then immediately re-listing in one invocation is rare enough (and cheap
        # enough as two commands) that it is not worth the ambiguity of picking an
        # implicit ordering between "mutate" and "report" in a single call.
        if ledger_flag:
            _usage_error("--mark-rebuilt and --ledger are mutually exclusive")
        if all_flag or dataset or regenerate or as_json:
            _usage_error(
                "--mark-rebuilt takes its own dataset name; it cannot be combined "
                "with --all, --regenerate, --json, or a dataset argument"
            )
        try:
            _mark_rebuilt(mark_rebuilt)
        except KeyError:
            click.echo(f"dataset not in ledger: {mark_rebuilt}", err=True)
            raise SystemExit(EXIT_REGISTRY_ERROR)
        except drift_ledger.LedgerError as exc:
            # Read leniently, a conflicted ledger looked empty and every name in it
            # "was not in the ledger"; writing would have erased it. Nothing is touched.
            click.echo(f"{exc}; nothing was written", err=True)
            raise SystemExit(EXIT_PROBE_FAILED)
        click.echo(f"marked rebuilt: {mark_rebuilt}")
        return

    if ledger_flag:
        # --ledger is its own reporting mode over the whole ledger, not a per-dataset
        # drift check -- --all / a dataset argument / --regenerate / --json all
        # describe a *different* command (run the probes; optionally on one dataset;
        # optionally overwriting the baseline; optionally as JSON). Silently accepting
        # them alongside --ledger and running the ledger dump anyway discards
        # whatever the other flag asked for without saying so; --domain is left
        # unguarded because it is already a no-op without --all outside this branch
        # too, so --ledger --domain is not a NEW inconsistency.
        if all_flag or dataset or regenerate or as_json:
            _usage_error(
                "--ledger cannot be combined with --all, --regenerate, --json, or a "
                "dataset argument"
            )
        # Strict, like --mark-rebuilt and --regenerate: a lenient read turned a merge
        # conflict or a wrong-shaped ledger into {}, so this used to report "no datasets
        # pending rebuild" for a file the writers already refuse to touch. A missing
        # file is still {} under strict=True, so a fresh checkout stays a clean sweep.
        try:
            ledger = drift_ledger.load_ledger(LEDGER_PATH, strict=True)
        except drift_ledger.LedgerError as exc:
            click.echo(str(exc), err=True)
            raise SystemExit(EXIT_PROBE_FAILED)
        stale = _stale_datasets(ledger)
        if not stale:
            click.echo("no datasets pending rebuild")
            return
        for name, entry in stale:
            click.echo(
                f"{name}\tupstream={entry['last_upstream_change']}\t"
                f"rebuilt={entry['rebuilt_at'] or 'never'}"
            )
        return

    if all_flag and dataset:
        _usage_error("pass either a dataset name or --all, not both")
    if not all_flag and not dataset:
        _usage_error("specify a dataset name or --all")

    reg = plugin_loader.get_registry()

    # A unit that failed to load contributes no rows at all -- not even a
    # `probe_failed` one -- so `--all` would report a clean sweep over a silently
    # smaller set. `drift-health.yml` greps the JSON for `status == "probe_failed"` and
    # gates on the exit code, so a stderr warning alone is read by humans and by nothing
    # else: the workflow still prints "all probes healthy" and opens no issue.
    #
    # It therefore becomes a real exit code below, and -- in `--json` mode, which is what
    # the workflows capture -- a real row too. The human-readable branch still gets only
    # the stderr warning above, which is the right split: a person reads the warning, CI
    # reads the row. Together they make the existing workflow file the issue it already
    # knows how to file (#351).
    load_errors = _relevant_load_errors(
        reg,
        reg.load_errors(),
        dataset=None if all_flag else dataset,
        domain=domain if all_flag else None,
    )
    for unit, exc in load_errors:
        click.echo(
            f"WARNING: {unit!r} failed to load, so it is NOT being drift-checked: {exc}",
            err=True,
        )

    if regenerate:
        if all_flag:
            _usage_error("--regenerate requires a specific dataset name")
        if as_json:
            _usage_error("--regenerate does not support --json")
        if not no_ledger_row:
            # Before the probe runs: the row is written after the fingerprint, and a
            # ledger that cannot be parsed must stop the regenerate while nothing has
            # been rewritten yet, not after.
            try:
                drift_ledger.load_ledger(LEDGER_PATH, strict=True)
            except drift_ledger.LedgerError as exc:
                click.echo(
                    f"{exc}; repair it before regenerating {dataset}. Nothing was "
                    "written.",
                    err=True,
                )
                raise SystemExit(EXIT_PROBE_FAILED)
        try:
            previous, observed = _regenerate_fingerprint(reg, dataset, timeout=timeout)
        except KeyError:
            click.echo(f"unknown dataset: {dataset}", err=True)
            _echo_unknown_dataset_hint(reg)
            raise SystemExit(EXIT_REGISTRY_ERROR)
        except DriftProbeError as exc:
            # The probe raised, timed out, or returned a non-mapping: nothing was written.
            # Exit 2, not a traceback's 1 (== EXIT_DRIFTED): drift_to_pr.py treats any
            # non-zero here as "do not commit", which is the right outcome either way, but
            # a human reading exit codes must not see "drifted" for "the probe broke".
            click.echo(
                f"probe failed for {dataset}; fingerprint NOT rewritten: {exc}",
                err=True,
            )
            raise SystemExit(EXIT_PROBE_FAILED)
        except OSError as exc:
            # The probe succeeded but write_fingerprint's write/replace failed (disk
            # full, permission denied, ...). Same outcome for the caller as a probe
            # failure -- nothing was rewritten -- so the same exit 2, not an unhandled
            # traceback's 1 (== EXIT_DRIFTED).
            click.echo(
                f"could not write the fingerprint for {dataset}: {exc}", err=True
            )
            raise SystemExit(EXIT_PROBE_FAILED)
        except PluginLoadError as exc:
            click.echo(
                f"{dataset} failed to load, so its fingerprint cannot be regenerated: "
                f"{exc}",
                err=True,
            )
            raise SystemExit(EXIT_PROBE_FAILED)
        except Exception as exc:  # noqa: BLE001
            # Anything else write_fingerprint or the probe can raise that is not one of
            # the specific cases above -- e.g. a probe returning a dict with a
            # non-string-safe key makes json.dumps raise TypeError deep inside
            # write_fingerprint. Before this, that escaped as an unhandled traceback's
            # exit 1 (== EXIT_DRIFTED), same failure mode #361 fixed for DriftProbeError
            # and OSError specifically.
            click.echo(
                f"could not regenerate the fingerprint for {dataset}: "
                f"{type(exc).__name__}: {exc}",
                err=True,
            )
            raise SystemExit(EXIT_PROBE_FAILED)
        click.echo(f"regenerated: {dataset}")
        if not no_ledger_row:
            # Before the load-error exit below: the new baseline is already on disk, so
            # its row belongs beside it whatever this run's exit code turns out to be.
            _record_ledger_row(
                dataset,
                previous,
                observed,
                accepted_in or drift_ledger.MANUAL_ACCEPTED_IN,
            )
        if load_errors:
            # The exit-code block below is never reached on this path, so the #351 guard
            # did not apply to --regenerate. Only units relevant to THIS dataset are in
            # `load_errors` now, so this fires when its provider is broken, not when an
            # unrelated plugin is (drift_to_pr.py runs this step unattended).
            click.echo(
                f"WARNING: {len(load_errors)} unit(s) related to {dataset} failed to "
                f"load (listed above); exiting {EXIT_PROBE_FAILED}",
                err=True,
            )
            raise SystemExit(EXIT_PROBE_FAILED)
        return

    try:
        targets = (
            reg.list_datasets(domain=domain) if all_flag else [reg.get_dataset(dataset)]
        )
    except KeyError:
        click.echo(f"unknown dataset: {dataset}", err=True)
        _echo_unknown_dataset_hint(reg)
        raise SystemExit(EXIT_REGISTRY_ERROR)
    except PluginLoadError as exc:
        # A dataset whose callables would not import: the registry cached the failure and
        # get_dataset re-raises it. Not a KeyError, so before #364 it escaped as a
        # traceback -- exit 1, which a wrapper reads as EXIT_DRIFTED.
        if as_json:
            click.echo(
                json.dumps([_load_error_row(dataset, exc)], indent=2, default=str)
            )
        click.echo(
            f"{dataset} failed to load, so it cannot be drift-checked: {exc}", err=True
        )
        raise SystemExit(EXIT_PROBE_FAILED)

    # An empty sweep is the limiting case of the silently-smaller set #351 is about, and
    # it slips past every guard above: no targets means no results, `exit_codes` stays
    # {EXIT_CLEAN}, and `drift.yml` (rc > 2) plus `drift_to_pr.py` both pass on a report
    # of `[]`. The workflow goes green having checked nothing.
    #
    # The two ways to get here need different answers. A --domain that no manifest
    # declares is a caller mistake -- the option takes a free-form string and its help
    # text names no valid value, so a typo silently means "check nothing" -- and a
    # UsageError can name the real domains. `known` comes from the manifests, bound or
    # not: a domain whose every dataset failed to bind is still a real domain, and
    # listing it as "(none)" misdescribed the registry. An empty registry with no filter
    # is an infrastructure failure (nothing discovered, nothing shipped), which is what
    # EXIT_PROBE_FAILED means.
    #
    # Neither applies when the sweep is empty BECAUSE things failed to load: the
    # manifests are registered, they did not bind, and `load_errors` holds a row for
    # each. Exiting here printed no JSON at all, so `drift-health.yml` `json.load`ed an
    # empty capture and went red on a parse error instead of filing the probe-failed
    # issue it knows how to file. Falling through emits those rows and exits
    # EXIT_PROBE_FAILED through the load-error rule below.
    if all_flag and not targets:
        known = sorted({m.domain for m in reg.list_manifests() if m.domain})
        if domain is not None and domain not in known:
            _usage_error(
                f"--domain {domain!r} matches no dataset; known domains: "
                f"{', '.join(known) or '(none)'}"
            )
        if not load_errors:
            click.echo(
                "no datasets registered, so nothing was drift-checked -- refusing to "
                "report a clean sweep over an empty set",
                err=True,
            )
            raise SystemExit(EXIT_PROBE_FAILED)

    # run_drift_checks, not a comprehension over run_drift_check: datasets sharing a
    # baseline AND a probe callable share one drift signal, so it is probed once and
    # fanned out. Every dataset still gets its own entry; they just agree, and the CI
    # bot collapses them into one PR via the shared fingerprint_path.
    results = drift_runner.run_drift_checks(targets, timeout=timeout)

    if as_json:
        rows = [_serialize(r) for r in results]
        rows.extend(_load_error_row(unit, exc) for unit, exc in load_errors)
        click.echo(json.dumps(rows, indent=2, default=str))
    else:
        for r in results:
            click.echo(f"{r.dataset_name}: {r.status}")
            if r.probe_error is not None:
                # The reason was serialised only under --json, so a person running
                # `hvantk drift x:y` saw `probe_failed` and nothing else. stderr, like
                # every other diagnostic here: stdout stays the per-dataset summary.
                click.echo(f"  {r.probe_error}", err=True)
            if r.diff:
                click.echo(json.dumps(r.diff, indent=2, default=str))

    # Emit stub WARNINGs to stderr in BOTH human-readable and --json modes.
    # The scheduled drift workflow captures `drift --all --json` stdout to a
    # file (.github/workflows/drift.yml), so a WARNING confined to the
    # human-readable path would never surface in CI logs — re-hiding doc-only
    # stubs. stderr keeps machine-readable stdout clean while CI logs stay
    # visibly non-green.
    for r in results:
        if r.status == "stub":
            reason = (r.observed or {}).get("reason", "no programmatic source")
            click.echo(
                f"WARNING: {r.dataset_name}: stub probe — {reason}; "
                "no real drift detection (documentation-only source).",
                err=True,
            )

    # Priority: probe_failed(2) > drifted(1) > clean(0). Infra failure trumps
    # drift because drifted output is only meaningful if the probe actually ran.
    # status="stub" is intentional (doc-only source) → exit clean, but the
    # WARNING above keeps it from being a silent false-green.
    exit_codes = {EXIT_CLEAN}
    for r in results:
        if r.status == "drifted":
            exit_codes.add(EXIT_DRIFTED)
        elif r.status == "probe_failed":
            exit_codes.add(EXIT_PROBE_FAILED)
    # A unit that never registered is a probe that cannot run. Reporting clean over a
    # silently smaller set is the failure mode #351 exists to stop, and `drift.yml` /
    # `drift-health.yml` both gate on `rc > 2`, so 2 still reads as "probe trouble"
    # rather than "the command itself broke".
    if load_errors:
        exit_codes.add(EXIT_PROBE_FAILED)
    raise SystemExit(max(exit_codes))


def _echo_unknown_dataset_hint(reg) -> None:
    """A second stderr line pointing at `hvantk plugins errors` when the registry has
    load failures that might be why the caller's dataset was not found.

    `unknown dataset: gwas-catalog:associations` alone reads as a typo, but the real
    cause can be that the whole `gwas_catalog` directory (dir name, not the
    hyphenated `name:` the manifest declares) never loaded at all -- the caller typed
    the provider's real name and got a KeyError with no mention of the load failure
    that explains it.
    """
    # A `cli:<command>` row is a command that did not import, never a missing dataset.
    n = sum(1 for unit, _ in reg.load_errors() if not plugin_loader.is_cli_unit(unit))
    if n:
        click.echo(
            f"({n} unit(s) failed to load and may include it; see: hvantk plugins errors)",
            err=True,
        )


def _regenerate_fingerprint(
    reg, dataset_name: str, *, timeout: int
) -> tuple[dict | None, dict]:
    """Regenerate ``dataset_name``'s baseline; return (the baseline it replaced, the
    new one). The old one is read first because the regenerate overwrites it, and the
    ledger row depends on what actually changed."""
    spec = reg.get_dataset(dataset_name)
    previous = drift_ledger.read_baseline(Path(spec.test_paths.drift_fingerprint))
    observed = drift_runner.regenerate_fingerprint(spec, timeout=timeout)
    return previous, observed


def _record_ledger_row(
    dataset_name: str, previous: dict | None, observed: dict, accepted_in: str
) -> None:
    """Record the rebuild-ledger row for a regenerate that captured an upstream change,
    or say why it recorded none (`drift_ledger.record_regeneration` decides).

    The fingerprint is already rewritten when this runs, so a failed ledger write must
    say so: a plain re-run would find the baseline unchanged and record nothing. Any
    failure exits EXIT_PROBE_FAILED, never a traceback's 1 (== EXIT_DRIFTED), for the
    reason the regenerate branch above gives.
    """
    try:
        row, skipped = drift_ledger.record_regeneration(
            dataset_name,
            previous=previous,
            regenerated=observed,
            accepted_in=accepted_in,
            path=LEDGER_PATH,
        )
    except Exception as exc:  # noqa: BLE001
        click.echo(
            f"the fingerprint for {dataset_name} was regenerated, but its rebuild-ledger "
            f"row could not be written to {LEDGER_PATH}: {exc}. Restore the old "
            "baseline and re-run, or the change goes unrecorded.",
            err=True,
        )
        raise SystemExit(EXIT_PROBE_FAILED) from exc
    if row is None:
        click.echo(f"no ledger row recorded for {dataset_name}: {skipped}")
        return
    click.echo(
        f"ledger row recorded for {dataset_name} (signal: {row['signal']}, "
        f"accepted_in: {row['accepted_in']}); commit {LEDGER_PATH} with the "
        "fingerprint"
    )


def _load_error_row(unit: str, exc: Exception) -> dict:
    """The JSON row a unit that never registered gets -- the same shape a failing probe
    produces, so no consumer needs to learn a new key."""
    return {
        "dataset_name": unit,
        "status": "probe_failed",
        "observed": None,
        "expected": None,
        "diff": None,
        "probe_error": f"plugin failed to load: {exc}",
        "fingerprint_path": None,
        "probe_ref": None,
    }


def _relevant_load_errors(
    reg, load_errors, *, dataset: str | None, domain: str | None
) -> list[tuple[str, Exception]]:
    """The load errors that bear on THIS invocation.

    Applying every recorded failure to every run made `hvantk drift clinvar:variants`
    report -- and fail on -- plugins the caller never asked about, and made every
    unattended `--regenerate` fail whenever any unrelated plugin was broken (#364).

    Units come in four shapes (`loader._record_load_error`): a bare provider name, a
    `provider:dataset` key, `entry-point:<name>`, or `cli:<command>`. A `cli:` unit is
    never relevant: it records a manifest command that did not import, and a dataset
    that fails with it has a row of its own (the reverted first fix, which filed it
    under the provider, made `drift` exit 2 with a bogus row in any process that had
    wired the downloaders, #408). For one dataset, its own key, its provider and its
    provider's entry point are relevant; for `--all --domain X`, a dataset-level unit
    whose manifest declares another domain is not, while provider and entry-point units
    carry no domain and are kept -- dropping them would recreate the silently-smaller
    sweep #351 exists to stop.

    A bare provider-level unit is matched against the requested dataset's provider after
    normalising ``_``/``-`` on both sides: `load_from_skills_root` records the plugin
    DIRECTORY's name when `plugin.yaml` is missing, and `_provider_id_hint` falls back to
    it too when the manifest's own `name` cannot be read -- and 9 of the 21 in-tree
    providers have a directory name that differs from the manifest's `name:` only by
    `_` vs `-` (`gwas_catalog` vs `gwas-catalog`, `uniprot_ptm` vs `uniprot-ptm`, ...).
    Exact string equality silently dropped those. The `entry-point:<name>` match stays an
    exact comparison -- it is keyed by the installed entry point's own name, which is
    expected to follow the provider's naming convention already.
    """
    load_errors = [
        (unit, exc) for unit, exc in load_errors if not plugin_loader.is_cli_unit(unit)
    ]
    if dataset is not None:
        provider = dataset.split(":", 1)[0]
        wanted = {dataset, f"entry-point:{provider}"}
        return [
            (unit, exc)
            for unit, exc in load_errors
            if unit in wanted or unit.replace("_", "-") == provider.replace("_", "-")
        ]
    if domain is None:
        return list(load_errors)
    domains = {m.name: m.domain for m in reg.list_manifests()}
    return [
        (unit, exc)
        for unit, exc in load_errors
        if unit.startswith("entry-point:")
        or ":" not in unit
        or domains.get(unit, domain) == domain
    ]


def _serialize(result: drift_runner.DriftResult) -> dict:
    d = dataclasses.asdict(result)
    if d["probe_error"] is not None:
        d["probe_error"] = str(result.probe_error)
    return d


def _load_ledger() -> dict:
    """Read the ledger at `LEDGER_PATH`. A missing or corrupt file -- or JSON that
    parses but is not an object -- yields {}, so `_stale_datasets` never calls
    `.items()` on a list (see `drift_ledger.load_ledger`)."""
    return drift_ledger.load_ledger(LEDGER_PATH)


def _parse_iso(s: str) -> datetime | None:
    """Parse an ISO-8601 timestamp (tolerating a trailing ``Z``) to a tz-aware
    ``datetime``. Returns None if ``s`` is not parseable. Mirrors the approach
    ``should_escalate`` uses in ``.github/scripts/drift_to_pr.py``.
    """
    try:
        return datetime.fromisoformat(s.replace("Z", "+00:00"))
    except (ValueError, AttributeError, TypeError):
        return None


def _stale_datasets(ledger: dict) -> list[tuple[str, dict]]:
    """Datasets whose upstream moved after their last rebuild. Never-rebuilt counts.

    Timestamps are parsed to real instants rather than compared as strings --
    "2026-08-23T09:00:00-05:00" (=14:00 UTC) sorts BEFORE
    "2026-08-23T10:00:00+00:00" lexicographically despite being the LATER instant,
    which would misreport a freshly-rebuilt dataset as stale. An unparseable
    timestamp on either side is treated as stale rather than silently as fresh: this
    is a "which artifacts are stale?" report, and a false "stale" merely prompts
    someone to look, while a false "fresh" would hide a genuinely stale artifact.
    """
    out = []
    for name, entry in sorted(ledger.items()):
        rebuilt = entry.get("rebuilt_at")
        if rebuilt is None:
            out.append((name, entry))
            continue
        rebuilt_at = _parse_iso(rebuilt)
        last_upstream_change = _parse_iso(entry.get("last_upstream_change", ""))
        if (
            rebuilt_at is None
            or last_upstream_change is None
            or rebuilt_at < last_upstream_change
        ):
            out.append((name, entry))
    return out


def _mark_rebuilt(dataset_name: str) -> None:
    """Set ``dataset_name``'s ``rebuilt_at`` to now in the ledger at `LEDGER_PATH`.

    Raises ``KeyError`` when the dataset has no row and `drift_ledger.LedgerError` when
    the ledger cannot be parsed; neither writes anything (see `drift_ledger.mark_rebuilt`).
    """
    drift_ledger.mark_rebuilt(dataset_name, path=LEDGER_PATH)
