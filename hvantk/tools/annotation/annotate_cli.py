"""CLI for building annotation artifacts.

    hvantk annotate spine --genes G.ht --structure S.ht --hgnc H.ht --output spine.ht

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
@click.option("--genes", required=True, help="Path to the ensembl-gene:genes table.")
@click.option(
    "--structure", required=True, help="Path to the ensembl-gene:structure table."
)
@click.option("--hgnc", required=True, help="Path to the hgnc:lookup table.")
@click.option("--output", required=True, help="Output path for the spine table.")
@click.option(
    "--biotype",
    default="protein_coding",
    show_default=True,
    help="Gene biotype to keep; pass 'all' to keep every biotype.",
)
def spine_cmd(genes, structure, hgnc, output, biotype):
    """Build the gene spine every annotation joins onto."""
    import hail as hl

    from hvantk.algorithms.annotation.spine import build_spine, spine_mapping_rate
    from hvantk.core.utils.hail_context import init_hail

    init_hail()

    ht = build_spine(
        hl.read_table(genes),
        hl.read_table(structure),
        hl.read_table(hgnc),
        biotype=None if biotype == "all" else biotype,
    )
    ht.write(output, overwrite=True)

    ht = hl.read_table(output)
    click.echo(f"spine: {ht.count()} genes -> {output}")
    click.echo(f"HGNC mapping rate: {spine_mapping_rate(ht):.2%}")
