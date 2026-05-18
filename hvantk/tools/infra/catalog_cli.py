"""Inspect the hvantk dataset catalog aggregated from per-plugin catalogs.

This CLI is a read-only window onto the live registry: it does NOT load the
old ``hvantk/resources/catalog.yaml`` summary file (removed when per-domain
``datasets.json`` files were retired). Instead it instantiates
:class:`hvantk.resources.unified_registry.HvantkRegistry`, which aggregates
each plugin's ``catalog/datasets.json`` plus the surviving legacy
``resources/registry/genomics/datasets.json``.

Subcommands:

* ``list``    — table view with omics/source/organism filters
* ``show``    — full record for a single accession (yaml or json)
* ``stats``   — aggregate counts by omics, organism, source
* ``search``  — free-text query over title / description / accession

The legacy ``catalog build`` subcommand was removed; use
``hvantk reprocess <provider:dataset>`` for the equivalent build workflow.
"""

from __future__ import annotations

import json
import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS
from hvantk.resources.unified_registry import HvantkRegistry

logger = logging.getLogger(__name__)

_OMICS_CHOICES = ["transcriptomics", "proteomics", "genomics", "epigenomics"]


@click.group(
    "catalog",
    context_settings=CONTEXT_SETTINGS,
    help=(
        "Inspect the hvantk dataset catalog (aggregated from per-plugin catalogs). "
        "Use `hvantk reprocess <provider:dataset>` to rebuild a dataset; "
        "the legacy `catalog build` subcommand was removed."
    ),
)
def catalog():
    """Root command group for catalog operations."""


@catalog.command("list")
@click.option(
    "--omics-type",
    type=click.Choice(_OMICS_CHOICES),
    default=None,
    help="Filter by omics type.",
)
@click.option(
    "--data-source",
    default=None,
    help="Filter by data_source field (substring match, case-insensitive).",
)
@click.option(
    "--organism",
    default=None,
    help="Filter by organism field (substring match, case-insensitive).",
)
@click.option(
    "--limit",
    default=50,
    show_default=True,
    type=int,
    help="Maximum number of rows to print.",
)
def list_entries(omics_type, data_source, organism, limit):
    """List catalog entries with optional filters."""
    reg = HvantkRegistry()
    entries = reg.search(
        query="",
        omics_type=omics_type,
        organism=organism,
        data_source=data_source,
    )

    if not entries:
        click.echo("(no entries match)")
        return

    click.echo(f"{'ACCESSION':<24} {'OMICS':<16} {'SOURCE':<20} TITLE")
    for entry in entries[:limit]:
        acc = entry.get("accession", "")
        omics = entry.get("_omics_type", "")
        src = entry.get("data_source", "")
        title = (entry.get("title") or "")[:60]
        click.echo(f"{acc:<24} {omics:<16} {src:<20} {title}")
    if len(entries) > limit:
        click.echo(f"... ({len(entries) - limit} more; raise --limit to see them)")


@catalog.command("show")
@click.argument("accession")
@click.option(
    "--format",
    "fmt",
    default="yaml",
    type=click.Choice(["yaml", "json"]),
    show_default=True,
    help="Output format.",
)
def show_entry(accession, fmt):
    """Show full record for one catalog entry by accession."""
    reg = HvantkRegistry()
    entry = reg.get_dataset(accession)
    if not entry:
        raise click.ClickException(f"Entry not found: {accession}")
    if fmt == "json":
        click.echo(json.dumps(entry, indent=2, default=str))
    else:
        import yaml

        click.echo(yaml.safe_dump(entry, sort_keys=False))


@catalog.command("stats")
def stats():
    """Aggregate catalog statistics (counts by omics, organism, source)."""
    reg = HvantkRegistry()
    s = reg.get_stats()

    click.echo(f"total datasets: {s['total_datasets']}")
    click.echo("")
    click.echo("by omics type:")
    for omics, count in s.get("by_omics_type", {}).items():
        click.echo(f"  {omics:<16} {count}")
    click.echo("")

    top_organisms = sorted(
        s.get("organisms", {}).items(), key=lambda x: -x[1]
    )[:10]
    click.echo("top organisms:")
    for org, count in top_organisms:
        click.echo(f"  {org:<40} {count}")
    click.echo("")

    top_sources = sorted(
        s.get("data_sources", {}).items(), key=lambda x: -x[1]
    )[:10]
    click.echo("top data sources:")
    for src, count in top_sources:
        click.echo(f"  {src:<40} {count}")


@catalog.command("search")
@click.argument("query")
@click.option(
    "--omics-type",
    type=click.Choice(_OMICS_CHOICES),
    default=None,
    help="Restrict the search to one omics bucket.",
)
@click.option(
    "--limit",
    default=50,
    show_default=True,
    type=int,
    help="Maximum number of rows to print.",
)
def search_entries(query, omics_type, limit):
    """Free-text search across title / description / accession."""
    reg = HvantkRegistry()
    entries = reg.search(query=query, omics_type=omics_type)
    if not entries:
        click.echo("(no matches)")
        return

    click.echo(f"matches: {len(entries)}")
    click.echo(f"{'ACCESSION':<24} {'OMICS':<16} TITLE")
    for entry in entries[:limit]:
        acc = entry.get("accession", "")
        omics = entry.get("_omics_type", "")
        title = (entry.get("title") or "")[:60]
        click.echo(f"{acc:<24} {omics:<16} {title}")
    if len(entries) > limit:
        click.echo(f"... ({len(entries) - limit} more; raise --limit)")
