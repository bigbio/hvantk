"""`hvantk plugins ...` commands for introspecting the plugin registry."""

from __future__ import annotations

import json
from pathlib import Path

import click

from hvantk.core.plugin import loader as plugin_loader


@click.group(name="plugins")
def plugins_group():
    """Inspect the hvantk plugin registry."""


@plugins_group.command(name="list")
def list_cmd():
    """List loaded provider plugins."""
    reg = plugin_loader.get_registry()
    providers = reg.list_providers()
    if not providers:
        click.echo("(no plugins loaded)")
        return
    click.echo(f"{'NAME':<24} {'VERSION':<12} {'DATASETS':<8}")
    for p in providers:
        click.echo(f"{p.name:<24} {p.version:<12} {len(p.datasets):<8}")


@plugins_group.command(name="describe")
@click.argument("provider")
def describe_cmd(provider: str):
    """Show full details for one provider."""
    reg = plugin_loader.get_registry()
    try:
        p = reg.get_provider(provider)
    except KeyError:
        raise click.ClickException(f"unknown provider: {provider}")
    click.echo(f"name:    {p.name}")
    click.echo(f"version: {p.version}")
    click.echo("datasets:")
    for ds in p.datasets:
        click.echo(f"  - {ds.name}  ({ds.domain}, {ds.backend})")
        click.echo(f"      skill: {ds.skill_path}")
        click.echo(f"      drift_fingerprint: {ds.test_paths.drift_fingerprint}")


@plugins_group.command(name="errors")
def errors_cmd():
    """List plugins that failed to load and why."""
    reg = plugin_loader.get_registry()
    errs = reg.load_errors()
    if not errs:
        click.echo("(no load errors)")
        return
    for plugin_id, exc in errs:
        click.echo(f"{plugin_id}: {exc}")


@plugins_group.command(name="validate")
@click.argument("manifest_path", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--strict-artifacts/--no-strict-artifacts",
    default=False,
    help=(
        "Fail if the manifest declares a fixture/snapshot/fingerprint that is not on "
        "disk. Warns by default, because most in-tree plugins do not yet ship a full "
        "set; see the KNOWN_INCOMPLETE ledger in the test suite."
    ),
)
def validate_cmd(manifest_path: str, strict_artifacts: bool):
    """Validate a plugin.yaml file against the schema (offline)."""
    import yaml
    import jsonschema

    schema = json.loads(
        (Path(plugin_loader.__file__).parent / "manifest.schema.json").read_text()
    )
    content = yaml.safe_load(Path(manifest_path).read_text())
    try:
        jsonschema.validate(content, schema)
    except jsonschema.ValidationError as exc:
        raise click.ClickException(f"validation failed: {exc.message}")
    catalog_rel = content.get("catalog")
    if catalog_rel:
        catalog_path = Path(manifest_path).parent / catalog_rel
        if not catalog_path.is_file():
            raise click.ClickException(f"catalog not found: {catalog_path}")
        try:
            catalog_schema = json.loads(
                (Path(plugin_loader.__file__).parent / "catalog_entry.schema.json").read_text()
            )
        except (OSError, json.JSONDecodeError) as exc:
            raise click.ClickException(
                f"failed to read bundled catalog schema: {exc}"
            ) from exc
        try:
            entries = json.loads(catalog_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise click.ClickException(
                f"failed to read catalog {catalog_path}: {exc}"
            ) from exc
        if not isinstance(entries, list):
            raise click.ClickException(f"catalog must be a JSON array: {catalog_path}")
        errors: list[str] = []
        seen: set[str] = set()
        for i, entry in enumerate(entries):
            try:
                jsonschema.validate(entry, catalog_schema)
            except jsonschema.ValidationError as exc:
                acc = entry.get("accession", f"index {i}") if isinstance(entry, dict) else f"index {i}"
                errors.append(f"entry {acc}: {exc.message}")
                continue
            acc = entry["accession"]
            if acc in seen:
                errors.append(f"duplicate accession within catalog: {acc}")
            seen.add(acc)
        if errors:
            raise click.ClickException(
                "catalog validation failed:\n  - " + "\n  - ".join(errors)
            )
        click.echo(f"catalog ok: {catalog_path} ({len(entries)} entries)")

    # Declared-but-absent validation artifacts. The loader resolves these paths without
    # requiring them, so a manifest can promise a snapshot it never shipped and still
    # load clean; surface that here instead of letting it pass silently.
    artifact_fields = ("fixture", "schema_snapshot", "row_snapshot", "drift_fingerprint")
    plugin_dir = Path(manifest_path).parent
    missing: list[str] = []
    for dataset in content.get("datasets", []):
        tests_block = dataset.get("tests") or {}
        for field in artifact_fields:
            rel = tests_block.get(field)
            if rel and not (plugin_dir / rel).exists():
                missing.append(f"{dataset.get('name', '?')}: {field} -> {rel}")
    if missing:
        if strict_artifacts:
            raise click.ClickException(
                "declared validation artifacts are missing:\n  - " + "\n  - ".join(missing)
            )
        click.echo(
            "warning: declared validation artifacts are missing "
            f"({len(missing)}; re-run with --strict-artifacts to fail):\n  - "
            + "\n  - ".join(missing),
            err=True,
        )

    # Datasets may deliberately SHARE a drift_fingerprint -- that is how a manifest says
    # "these have one drift signal, not N", and the drift runner groups on it so one
    # upstream event yields one PR rather than one per dataset. Sharing is only coherent
    # when the probe is also shared: two datasets pointing different probes at one
    # baseline would each overwrite the other's, and whichever regenerated last would
    # define "clean" for both. That is a manifest error, so name it.
    by_fingerprint: dict[str, list[tuple[str, tuple[str, str]]]] = {}
    for dataset in content.get("datasets", []):
        rel = (dataset.get("tests") or {}).get("drift_fingerprint")
        probe = dataset.get("drift_probe") or {}
        if not rel:
            continue
        by_fingerprint.setdefault(rel, []).append(
            (dataset.get("name", "?"), (probe.get("module", ""), probe.get("function", "")))
        )
    conflicts = [
        f"{rel} is shared by "
        + ", ".join(f"{name} -> {mod}:{fn}" for name, (mod, fn) in members)
        for rel, members in sorted(by_fingerprint.items())
        if len({probe for _, probe in members}) > 1
    ]
    if conflicts:
        raise click.ClickException(
            "datasets share a drift_fingerprint but declare different drift probes; "
            "they would overwrite each other's baseline:\n  - " + "\n  - ".join(conflicts)
        )

    click.echo(f"ok: {manifest_path}")
