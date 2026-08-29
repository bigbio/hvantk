"""`hvantk drift ...` command for fingerprint-based drift detection."""

from __future__ import annotations

import dataclasses
import json
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
def drift_cmd(dataset, all_flag, domain, as_json, regenerate, timeout, ledger_flag):
    """Compare a plugin's live drift-probe fingerprint against the expected file."""
    if ledger_flag:
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


def _stale_datasets(ledger: dict) -> list[tuple[str, dict]]:
    """Datasets whose upstream moved after their last rebuild. Never-rebuilt counts."""
    out = []
    for name, entry in sorted(ledger.items()):
        rebuilt = entry.get("rebuilt_at")
        if rebuilt is None or rebuilt < entry.get("last_upstream_change", ""):
            out.append((name, entry))
    return out
