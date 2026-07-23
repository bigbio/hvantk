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
@click.option(
    "--hgnc",
    "hgnc_path",
    default=None,
    help="Path to the hgnc:lookup table (.ht); required when the entry's key is hgnc_id or symbol.",
)
def prepare_cmd(spec, axis, input_path, spine, output, hgnc_path):
    """Map one source onto the spine's gene_id and select its declared columns."""
    import hail as hl

    from hvantk.algorithms.annotation.mapping import enforce_rate
    from hvantk.algorithms.annotation.prepare import (
        prepare_matrix_source,
        prepare_source,
        prepare_variant_source,
    )
    from hvantk.algorithms.annotation.spec import load_spec
    from hvantk.core.utils.hail_context import init_hail

    init_hail()

    entry = load_spec(spec).entry(axis)
    spine_ids = hl.read_table(spine).gene_id.collect()

    needs_hgnc = entry.key in ("hgnc_id", "symbol") or (
        entry.key == "variant"
        and entry.aggregate is not None
        and entry.aggregate.to != "gene_id"
    )
    hgnc = None
    if needs_hgnc:
        if not hgnc_path:
            raise click.UsageError(
                f"entry {entry.source!r} needs --hgnc <hgnc:lookup .ht>"
            )
        from hvantk.skills.hgnc.streamers import HGNCGeneCatalogStreamer

        hgnc = HGNCGeneCatalogStreamer.from_path(hgnc_path)

    if entry.matrix is not None:
        import anndata

        matrix_ad = anndata.read_h5ad(input_path)
        prepared, report = prepare_matrix_source(matrix_ad, spine_ids, entry, hgnc=hgnc)
    elif entry.key == "variant":
        source_ht = hl.read_table(input_path)
        prepared, report = prepare_variant_source(
            source_ht, spine_ids, entry, hgnc=hgnc
        )
    else:
        source_ht = hl.read_table(input_path)
        prepared, report = prepare_source(source_ht, spine_ids, entry, hgnc=hgnc)
    enforce_rate(report, entry.min_mapping_rate)
    prepared.write(output, overwrite=True)
    click.echo(f"prepare {axis} ({entry.source}): {report.summary()} -> {output}")
