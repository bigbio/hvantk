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
def validate_cmd(manifest_path: str):
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
    click.echo(f"ok: {manifest_path}")
