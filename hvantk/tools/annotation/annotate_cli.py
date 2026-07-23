"""CLI for building annotation artifacts.

    hvantk annotate spine --gene-table G.ht --hgnc H.ht --output spine.ht

Later phases add ``prepare``, ``compose`` and ``cohort`` to this group.

The command does its work in-process. Scheduling is deliberately external -- sbatch
wrappers submit it -- so that hvantk never imports a scheduler.
"""
from __future__ import annotations

import logging

import click

logger = logging.getLogger(__name__)


@click.group("annotate")
def annotate_group():
    """Build gene- and cohort-level annotation artifacts."""


@annotate_group.command("spine")
@click.option(
    "--gene-table",
    "gene_table",
    required=True,
    help="Path to the ensembl-gene:structure table (the canonical gene table).",
)
@click.option("--hgnc", required=True, help="Path to the hgnc:lookup table.")
@click.option("--output", required=True, help="Output path for the spine table.")
@click.option(
    "--biotype",
    default="protein_coding",
    show_default=True,
    help="Gene biotype to keep; pass 'all' to keep every biotype.",
)
def spine_cmd(gene_table, hgnc, output, biotype):
    """Build the gene spine every annotation joins onto."""
    import hail as hl

    from hvantk.algorithms.annotation.spine import build_spine, spine_mapping_rate
    from hvantk.core.utils.hail_context import init_hail

    init_hail()

    ht = build_spine(
        hl.read_table(gene_table),
        hl.read_table(hgnc),
        biotype=None if biotype == "all" else biotype,
    )
    ht.write(output, overwrite=True)

    ht = hl.read_table(output)
    click.echo(f"spine: {ht.count()} genes -> {output}")
    click.echo(f"HGNC mapping rate: {spine_mapping_rate(ht):.2%}")


@annotate_group.command("prepare")
@click.option("--spec", required=True, help="Path to the feature-spec YAML.")
@click.option("--axis", required=True, help="Which layer1 axis entry to prepare.")
@click.option(
    "--input",
    "input_path",
    required=True,
    help="Path to the built source AnnotationTable (.ht).",
)
@click.option("--spine", required=True, help="Path to the gene spine table (.ht).")
@click.option(
    "--output", required=True, help="Output path for the prepared table (.ht)."
)
def prepare_cmd(spec, axis, input_path, spine, output):
    """Map one source onto the spine's gene_id and select its declared columns."""
    import hail as hl

    from hvantk.algorithms.annotation.mapping import enforce_rate
    from hvantk.algorithms.annotation.prepare import prepare_source
    from hvantk.algorithms.annotation.spec import load_spec
    from hvantk.core.utils.hail_context import init_hail

    init_hail()

    entry = load_spec(spec).entry(axis)
    source_ht = hl.read_table(input_path)
    spine_ids = hl.read_table(spine).gene_id.collect()

    prepared, report = prepare_source(source_ht, spine_ids, entry)
    enforce_rate(report, entry.min_mapping_rate)  # raises MappingRateError if too low
    prepared.write(output, overwrite=True)

    click.echo(f"prepare {axis} ({entry.source}): {report.summary()} -> {output}")
