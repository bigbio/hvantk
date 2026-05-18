"""`hvantk tools ...` commands for introspecting the tool manifest registry."""

from __future__ import annotations

import json
from pathlib import Path

import click

from hvantk.core import tool_loader


@click.group(name="tools")
def tools_group():
    """Inspect the hvantk tool registry."""


@tools_group.command(name="list")
@click.option("--domain", default=None, help="Filter by domain")
def list_cmd(domain):
    """List registered tools (those with a tool.yaml manifest)."""
    reg = tool_loader.get_registry()
    tools = reg.list_tools(domain=domain)
    if not tools:
        click.echo("(no tools registered)")
        return
    click.echo(f"{'NAME':<24} {'DOMAIN':<12} {'TYPE':<16} DESCRIPTION")
    for t in tools:
        click.echo(f"{t.name:<24} {t.domain:<12} {t.type:<16} {t.description}")


@tools_group.command(name="describe")
@click.argument("name")
def describe_cmd(name):
    """Show full details for one tool."""
    reg = tool_loader.get_registry()
    try:
        t = reg.get_tool(name)
    except KeyError:
        raise click.ClickException(f"unknown tool: {name}")
    click.echo(f"name:        {t.name}")
    click.echo(f"domain:      {t.domain}")
    click.echo(f"type:        {t.type}")
    click.echo(f"description: {t.description}")
    click.echo(f"cli:         {t.cli_module}.{t.cli_callable}")
    click.echo(f"requires:    hail={t.requires.hail}, network={t.requires.network}")
    click.echo(f"purpose:     {t.purpose_short}")
    if t.purpose_long:
        click.echo("")
        click.echo(t.purpose_long.rstrip())
    if t.subcommands:
        click.echo("")
        click.echo("subcommands:")
        for s in t.subcommands:
            click.echo(f"  - {s.name}: {s.description}")


@tools_group.command(name="errors")
def errors_cmd():
    """List tool manifests that failed to load and why."""
    reg = tool_loader.get_registry()
    errs = reg.load_errors()
    if not errs:
        click.echo("(no load errors)")
        return
    for manifest_path, exc in errs:
        click.echo(f"{manifest_path}: {exc}")


@tools_group.command(name="validate")
@click.argument("manifest_path", type=click.Path(exists=True, dir_okay=False))
def validate_cmd(manifest_path):
    """Validate one tool.yaml against the schema (offline)."""
    import yaml
    import jsonschema

    schema = json.loads(
        (Path(tool_loader.__file__).parent / "tool_manifest.schema.json").read_text()
    )
    content = yaml.safe_load(Path(manifest_path).read_text())
    try:
        jsonschema.validate(content, schema)
    except jsonschema.ValidationError as exc:
        raise click.ClickException(f"validation failed: {exc.message}")
    click.echo(f"ok: {manifest_path}")
