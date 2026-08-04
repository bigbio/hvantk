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


@click.command(name="drift")
@click.argument("dataset", required=False)
@click.option("--all", "all_flag", is_flag=True, help="Run drift check for every dataset")
@click.option("--domain", default=None, help="Filter by domain (with --all)")
@click.option("--json", "as_json", is_flag=True, help="Emit machine-readable JSON")
@click.option("--regenerate", is_flag=True, help="Overwrite drift_fingerprint.json with observed probe output")
@click.option("--timeout", default=60, show_default=True, help="Probe timeout in seconds")
def drift_cmd(dataset, all_flag, domain, as_json, regenerate, timeout):
    """Compare a plugin's live drift-probe fingerprint against the expected file."""
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
