"""`hvantk drift ...` command for fingerprint-based drift detection."""

from __future__ import annotations

import dataclasses
import json
from datetime import datetime, timezone
from pathlib import Path

import click

from hvantk.core.plugin import drift_runner, loader as plugin_loader


EXIT_CLEAN = 0
EXIT_DRIFTED = 1
EXIT_PROBE_FAILED = 2
EXIT_REGISTRY_ERROR = 3

# hvantk/tools/plugins/drift_cli.py -> parents[2] is hvantk/, so this resolves to
# hvantk/resources/drift_ledger.json -- the ledger the drift bot writes in the same
# commit as the fingerprints it describes (see .github/scripts/drift_to_pr.py).
LEDGER_PATH = Path(__file__).resolve().parents[2] / "resources" / "drift_ledger.json"


@click.command(name="drift")
@click.argument("dataset", required=False)
@click.option("--all", "all_flag", is_flag=True, help="Run drift check for every dataset")
@click.option("--domain", default=None, help="Filter by domain (with --all)")
@click.option("--json", "as_json", is_flag=True, help="Emit machine-readable JSON")
@click.option("--regenerate", is_flag=True, help="Overwrite drift_fingerprint.json with observed probe output")
@click.option("--timeout", default=60, show_default=True, help="Probe timeout in seconds")
@click.option("--ledger", "ledger_flag", is_flag=True,
              help="List datasets whose upstream moved since their last rebuild")
@click.option("--mark-rebuilt", "mark_rebuilt", default=None, metavar="DATASET",
              help="Record DATASET as rebuilt now, clearing it from --ledger's stale list")
def drift_cmd(dataset, all_flag, domain, as_json, regenerate, timeout, ledger_flag, mark_rebuilt):
    """Compare a plugin's live drift-probe fingerprint against the expected file."""
    if mark_rebuilt is not None:
        # A standalone action, like --regenerate: it takes its OWN dataset name, so a
        # positional `dataset` / --all / --regenerate alongside it would be ambiguous
        # about which dataset is meant. Mutually exclusive with --ledger too -- marking
        # and then immediately re-listing in one invocation is rare enough (and cheap
        # enough as two commands) that it is not worth the ambiguity of picking an
        # implicit ordering between "mutate" and "report" in a single call.
        if ledger_flag:
            raise click.UsageError("--mark-rebuilt and --ledger are mutually exclusive")
        if all_flag or dataset or regenerate or as_json:
            raise click.UsageError(
                "--mark-rebuilt takes its own dataset name; it cannot be combined "
                "with --all, --regenerate, --json, or a dataset argument"
            )
        try:
            _mark_rebuilt(mark_rebuilt)
        except KeyError:
            click.echo(f"dataset not in ledger: {mark_rebuilt}", err=True)
            raise SystemExit(EXIT_REGISTRY_ERROR)
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
            raise click.UsageError(
                "--ledger cannot be combined with --all, --regenerate, --json, or a "
                "dataset argument"
            )
        stale = _stale_datasets(_load_ledger())
        if not stale:
            click.echo("no datasets pending rebuild")
            return
        for name, entry in stale:
            click.echo(f"{name}\tupstream={entry['last_upstream_change']}\t"
                       f"rebuilt={entry['rebuilt_at'] or 'never'}")
        return

    if all_flag and dataset:
        raise click.UsageError("pass either a dataset name or --all, not both")
    if not all_flag and not dataset:
        raise click.UsageError("specify a dataset name or --all")

    reg = plugin_loader.get_registry()

    if regenerate:
        if all_flag:
            raise click.UsageError("--regenerate requires a specific dataset name")
        try:
            _regenerate_fingerprint(reg, dataset)
        except KeyError:
            click.echo(f"unknown dataset: {dataset}", err=True)
            raise SystemExit(EXIT_REGISTRY_ERROR)
        click.echo(f"regenerated: {dataset}")
        return

    try:
        targets = (
            reg.list_datasets(domain=domain) if all_flag else [reg.get_dataset(dataset)]
        )
    except KeyError:
        click.echo(f"unknown dataset: {dataset}", err=True)
        raise SystemExit(EXIT_REGISTRY_ERROR)

    # run_drift_checks, not a comprehension over run_drift_check: datasets sharing a
    # baseline AND a probe callable share one drift signal, so it is probed once and
    # fanned out. Every dataset still gets its own entry; they just agree, and the CI
    # bot collapses them into one PR via the shared fingerprint_path.
    results = drift_runner.run_drift_checks(targets, timeout=timeout)

    if as_json:
        click.echo(json.dumps([_serialize(r) for r in results], indent=2, default=str))
    else:
        for r in results:
            click.echo(f"{r.dataset_name}: {r.status}")
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
    raise SystemExit(max(exit_codes))


def _regenerate_fingerprint(reg, dataset_name: str) -> None:
    spec = reg.get_dataset(dataset_name)
    observed = spec.drift_probe()
    Path(spec.test_paths.drift_fingerprint).write_text(
        json.dumps(observed, indent=2, default=str)
    )


def _serialize(result: drift_runner.DriftResult) -> dict:
    d = dataclasses.asdict(result)
    if d["probe_error"] is not None:
        d["probe_error"] = str(result.probe_error)
    return d


def _load_ledger() -> dict:
    """Read the ledger. A missing or corrupt file -- or JSON that parses but is not an
    object (a truthy list, a bare string) -- yields {} rather than raising.
    `json.loads(text) or {}` looks like it covers this, but `or` only substitutes on a
    FALSY parse (`[]`, `0`, `""`, `null`); a populated list or non-empty string is
    truthy and would pass straight through to `_stale_datasets`, which calls
    `.items()` and raises `AttributeError` on anything that isn't a dict.
    """
    if not LEDGER_PATH.is_file():
        return {}
    try:
        data = json.loads(LEDGER_PATH.read_text())
    except ValueError:
        return {}
    return data if isinstance(data, dict) else {}


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
        if rebuilt_at is None or last_upstream_change is None or rebuilt_at < last_upstream_change:
            out.append((name, entry))
    return out


def _mark_rebuilt(dataset_name: str) -> None:
    """Set ``dataset_name``'s ledger entry ``rebuilt_at`` to now (UTC, ISO-8601) and
    write the ledger back. Every other entry is copied through untouched.

    Raises ``KeyError`` if ``dataset_name`` has no ledger row -- the ledger is written
    exclusively by the drift bot (``.github/scripts/drift_to_pr.py``) when a dataset
    actually drifts, so a name absent from it never had a pending-rebuild signal to
    clear. Silently creating a row here would let `--ledger` under-report just as
    easily as the string-comparison and truthy-JSON bugs this same file fixes elsewhere
    -- fail loudly instead.
    """
    ledger = _load_ledger()
    if dataset_name not in ledger:
        raise KeyError(dataset_name)
    updated = dict(ledger)
    updated[dataset_name] = {
        **ledger[dataset_name],
        "rebuilt_at": datetime.now(timezone.utc).isoformat(),
    }
    LEDGER_PATH.write_text(json.dumps(updated, indent=2, sort_keys=True) + "\n")
