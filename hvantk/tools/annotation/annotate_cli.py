"""CLI for building annotation artifacts.

    hvantk annotate spine --gene-table G.ht --hgnc H.ht --output spine.ht

``prepare`` and ``compose`` live in this group; ``cohort`` is a deliberately separate
top-level group (``hvantk cohort``) rather than a subcommand here.

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
    help="Path to the built source: an AnnotationTable (.ht), or an AnnData (.h5ad) "
    "for a matrix entry.",
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


def _parse_prepared_option(values: tuple[str, ...]) -> dict[str, str]:
    """Parse repeatable ``--prepared axis=path`` options into a dict.

    Pure Python -- no Hail, no spec -- so it is unit-testable without a Hail session.
    """
    prepared_paths: dict[str, str] = {}
    for value in values:
        axis, sep, path = value.partition("=")
        axis, path = axis.strip(), path.strip()
        if not sep or not axis or not path:
            raise click.UsageError(
                f"--prepared value {value!r} is malformed; expected axis=path"
            )
        if axis in prepared_paths:
            raise click.UsageError(
                f"--prepared given twice for axis {axis!r} "
                f"({prepared_paths[axis]!r} then {path!r}); pass each axis once"
            )
        prepared_paths[axis] = path
    return prepared_paths


def _validate_prepared_axes(spec, prepared_paths: dict[str, str]) -> None:
    """Raise a clear ``click.UsageError`` if ``--prepared`` doesn't match the spec's axes.

    Checked in pure Python, before ``compose`` (which raises a raw ``KeyError`` on a
    missing axis) and before any Hail call, so this is unit-testable without Hail.
    """
    spec_axes = {entry.axis for entry in spec.layer1}
    prepared_axes = set(prepared_paths)

    missing = sorted(spec_axes - prepared_axes)
    if missing:
        raise click.UsageError(
            "missing --prepared for spec axis(es): " + ", ".join(missing)
        )

    unknown = sorted(prepared_axes - spec_axes)
    if unknown:
        raise click.UsageError(
            "--prepared names axis(es) not declared in the spec "
            f"{spec.name!r}: " + ", ".join(unknown)
        )


@annotate_group.command("compose")
@click.option("--spec", required=True, help="Path to the feature-spec YAML.")
@click.option("--spine", required=True, help="Path to the gene spine table (.ht).")
@click.option(
    "--prepared",
    "prepared_opts",
    required=True,
    multiple=True,
    help="Repeatable axis=path, one per spec layer1 axis (e.g. constraint=constraint.ht).",
)
@click.option(
    "--output", required=True, help="Output path for the composed matrix (.ht)."
)
@click.option(
    "--manifest",
    default=None,
    help="Output path for the manifest JSON (default: <output>.manifest.json).",
)
def compose_cmd(spec, spine, prepared_opts, output, manifest):
    """Left-join every spec axis's prepared table onto the spine."""
    import json

    from hvantk.algorithms.annotation.compose import compose
    from hvantk.algorithms.annotation.spec import load_spec

    feature_spec = load_spec(spec)
    prepared_paths = _parse_prepared_option(prepared_opts)
    _validate_prepared_axes(feature_spec, prepared_paths)

    manifest_path = manifest or f"{output}.manifest.json"

    import hail as hl

    from hvantk.core.utils.hail_context import init_hail

    init_hail()

    spine_ht = hl.read_table(spine)
    prepared_by_axis = {
        axis: hl.read_table(path) for axis, path in prepared_paths.items()
    }

    composed, report = compose(spine_ht, prepared_by_axis, feature_spec)
    composed.write(output, overwrite=True)

    with open(manifest_path, "w") as fh:
        json.dump(report, fh, indent=2)

    click.echo(
        f"compose {feature_spec.name}: {report['n_genes']} genes x "
        f"{len(feature_spec.layer1)} axes -> {output} (manifest: {manifest_path})"
    )
