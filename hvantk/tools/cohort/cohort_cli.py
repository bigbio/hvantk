"""CLI for presenting an external cohort to hvantk.

    hvantk cohort validate --cohort cohort.yaml
    hvantk cohort attach   --cohort cohort.yaml --layer1 composed.ht --output attached.ht

Every cohort verb lives in this one group. A later increment adds `cohort reduce`, which
derives a gene-level cohort table from variant-level genotypes -- it produces exactly the
table this contract already describes, so nothing downstream changes.

Validation is pure Python and runs before any Hail call, so a malformed manifest fails
fast with a clear message instead of after a Spark start-up.
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path

import click

logger = logging.getLogger(__name__)


@click.group("cohort")
def cohort_group():
    """Validate and attach external cohorts."""


def _read_header(table_path: str) -> list[str]:
    """Return the column names of a delimited cohort table.

    Pure Python and header-only, so `validate` stays cheap and Hail-free.
    """
    path = Path(table_path)
    if not path.exists():
        raise click.UsageError(f"cohort table not found: {table_path}")

    delimiter = "\t" if path.suffix in (".tsv", ".txt") else ","
    with path.open() as fh:
        try:
            header = next(csv.reader(fh, delimiter=delimiter))
        except StopIteration:
            raise click.UsageError(f"cohort table is empty: {table_path}")
    return [h.strip() for h in header]


def _check_declared_columns_exist(manifest, header: list[str]) -> None:
    """Every declared column must exist in the table.

    A missing declared column must never become a silently absent feature -- that is
    how a mapping bug turns into unexplained model degradation.
    """
    missing = [c for c in manifest.declared_columns() if c not in header]
    if missing:
        raise click.UsageError(
            f"cohort table {manifest.table} is missing declared column(s): "
            + ", ".join(missing)
        )


def _load_manifest(cohort_path: str):
    """Parse a manifest, converting contract violations into clean CLI errors."""
    from hvantk.algorithms.cohort.spec import load_cohort

    try:
        return load_cohort(cohort_path)
    except click.UsageError:
        raise
    except Exception as exc:
        raise click.UsageError(f"invalid cohort manifest {cohort_path}: {exc}")


@cohort_group.command("validate")
@click.option("--cohort", required=True, help="Path to the cohort manifest YAML.")
def validate_cmd(cohort):
    """Check a cohort manifest against its table (and labels) without building anything."""
    manifest = _load_manifest(cohort)
    header = _read_header(manifest.table)
    _check_declared_columns_exist(manifest, header)

    click.echo(f"cohort {manifest.name}: key={manifest.key} table={manifest.table}")
    click.echo(f"  prior: {manifest.prior.column} ({manifest.prior.direction})")
    for entry in manifest.cohort_axes:
        click.echo(f"  axis {entry.axis}: {', '.join(entry.columns)}")

    if manifest.labels is not None:
        from hvantk.algorithms.cohort.labels import load_label_genes

        try:
            genes = load_label_genes(manifest.labels)
        except Exception as exc:
            raise click.UsageError(f"labels: {exc}")
        click.echo(f"  labels: {len(genes)} genes from {manifest.labels.gene_set}")
    else:
        click.echo("  labels: none declared (rerank will require them)")

    click.echo(f"OK: {len(manifest.declared_columns())} declared columns present")


@cohort_group.command("attach")
@click.option("--cohort", required=True, help="Path to the cohort manifest YAML.")
@click.option(
    "--layer1", required=True, help="Path to the composed Layer-1 matrix (.ht)."
)
@click.option(
    "--output", required=True, help="Output path for the attached table (.ht)."
)
@click.option(
    "--report",
    default=None,
    help="Output path for the report JSON (default: <output>.report.json).",
)
@click.option(
    "--hgnc",
    "hgnc_path",
    default=None,
    help="Path to the hgnc:lookup table (.ht); required when the cohort key is "
    "hgnc_id or symbol.",
)
def attach_cmd(cohort, layer1, output, report, hgnc_path):
    """Join a cohort's declared columns onto the Layer-1 matrix over the full spine."""
    import json

    manifest = _load_manifest(cohort)
    header = _read_header(manifest.table)
    _check_declared_columns_exist(manifest, header)

    # Labels are symbol-keyed regardless of the cohort's own key, so declaring labels
    # needs the catalog even for a gene_id-keyed cohort.
    needs_hgnc = manifest.key in ("hgnc_id", "symbol") or manifest.labels is not None
    if needs_hgnc and not hgnc_path:
        reason = (
            f"cohort key {manifest.key!r}"
            if manifest.key in ("hgnc_id", "symbol")
            else "declared labels (label gene sets are symbol-keyed)"
        )
        raise click.UsageError(f"{reason} needs --hgnc <hgnc:lookup .ht>")

    report_path = report or f"{output}.report.json"

    import hail as hl

    from hvantk.algorithms.cohort.attach import attach
    from hvantk.core.utils.hail_context import init_hail

    init_hail()

    hgnc = None
    if needs_hgnc:
        from hvantk.skills.hgnc.streamers import HGNCGeneCatalogStreamer

        hgnc = HGNCGeneCatalogStreamer.from_path(hgnc_path)

    delimiter = "\t" if Path(manifest.table).suffix in (".tsv", ".txt") else ","
    cohort_ht = hl.import_table(
        manifest.table, delimiter=delimiter, impute=True, key=manifest.key
    )

    attached, result = attach(hl.read_table(layer1), cohort_ht, manifest, hgnc=hgnc)
    attached.write(output, overwrite=True)

    with open(report_path, "w") as fh:
        json.dump(result, fh, indent=2)

    click.echo(
        f"attach {manifest.name}: {result['n_tested']}/{result['n_genes']} genes "
        f"tested -> {output} (report: {report_path})"
    )
