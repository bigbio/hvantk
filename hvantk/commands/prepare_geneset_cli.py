"""CLI command for preparing gene set collections from plain TSV files.

Converts a headerless two-column TSV (gene_set_name<TAB>gene_symbol) into
a GeneSetCollection JSON file, with optional HGNC symbol validation.
"""

import logging
import os
import sys
from pathlib import Path

import click

logger = logging.getLogger(__name__)


@click.command(name="prepare-geneset")
@click.option(
    "-i",
    "--input",
    "input_path",
    type=click.Path(exists=True),
    required=True,
    help=(
        "Input TSV file (headerless, two columns: " "gene_set_name<TAB>gene_symbol)."
    ),
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    required=True,
    help="Output path for GeneSetCollection JSON.",
)
@click.option(
    "--hgnc",
    type=click.Path(exists=True),
    default=None,
    help=(
        "Path to HGNC data (.ht or .tsv) for symbol validation " "and alias resolution."
    ),
)
@click.option(
    "--min-genes",
    type=click.IntRange(min=0),
    default=0,
    help="Exclude gene sets with fewer than this many genes. [default: 0]",
)
@click.option(
    "--background",
    type=click.Path(exists=True),
    default=None,
    help=(
        "Text file with background gene universe (one symbol per line). "
        "Default: union of all genes across sets."
    ),
)
@click.option(
    "--export-gmt",
    type=click.Path(),
    default=None,
    help="Additionally export as GMT file.",
)
@click.option("--overwrite", is_flag=True, help="Overwrite output if exists.")
def prepare_geneset_cmd(
    input_path,
    output,
    hgnc,
    min_genes,
    background,
    export_gmt,
    overwrite,
):
    """Prepare a GeneSetCollection JSON from a plain TSV file.

    Reads a headerless two-column TSV (gene_set_name<TAB>gene_symbol) and
    produces a GeneSetCollection JSON file compatible with:

    \b
      hvantk psroc --gene-sets OUTPUT
      hvantk enrichex burden -s OUTPUT
      hvantk enrichex overlap -s OUTPUT

    Optionally validates symbols against HGNC, resolving aliases and
    reporting unrecognized entries.

    \b
    Example:
      hvantk prepare-geneset -i panels.tsv -o panels.json --hgnc /data/hgnc.ht
    """
    from hvantk.utils.gene_sets import load_gene_set, load_gene_sets_from_dict
    from hvantk.utils.geneset_io import parse_geneset_tsv, validate_with_hgnc

    output_path = Path(output)
    if output_path.exists() and not overwrite:
        click.echo(
            f"Error: output file '{output}' already exists. "
            "Use --overwrite to replace.",
            err=True,
        )
        sys.exit(1)

    # 1. Parse input TSV.
    click.echo(f"Loading gene sets from {input_path} ...", err=True)
    try:
        result = parse_geneset_tsv(Path(input_path))
    except (ValueError, FileNotFoundError) as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    gene_sets = result.gene_sets
    for w in result.warnings:
        click.echo(f"  Warning: {w}", err=True)

    n_genes = len({g for genes in gene_sets.values() for g in genes})
    click.echo(
        f"Parsed {n_genes} unique genes across {len(gene_sets)} gene sets.",
        err=True,
    )

    # 2. Optional HGNC validation.
    hgnc_validated = False
    aliases_resolved = {}
    unrecognized_symbols = []

    if hgnc:
        click.echo(f"\nValidating against HGNC ({hgnc}) ...", err=True)
        vr = validate_with_hgnc(gene_sets, hgnc)
        gene_sets = vr.gene_sets
        hgnc_validated = True
        aliases_resolved = vr.aliases_resolved
        unrecognized_symbols = sorted(vr.unrecognized)

        # Per-set summary.
        for set_name, genes in gene_sets.items():
            n_set = len(genes)
            set_aliases = {
                old: new for old, new in aliases_resolved.items() if new in set(genes)
            }
            n_aliases = len(set_aliases)
            set_unrec = [g for g in genes if g in vr.unrecognized]

            line = f"  {set_name}: {n_set} genes"
            if n_aliases:
                mappings = ", ".join(
                    f"{old} -> {new}" for old, new in set_aliases.items()
                )
                line += f", {n_aliases} aliases resolved ({mappings})"
            if set_unrec:
                line += f", {len(set_unrec)} unrecognized"
            click.echo(line, err=True)

        # Overall summary.
        n_total = len({g for genes in gene_sets.values() for g in genes})
        n_recognized = len(vr.recognized)
        pct = (n_recognized / n_total * 100) if n_total else 0
        click.echo(
            f"\nSummary:\n"
            f"  Gene sets:          {len(gene_sets)}\n"
            f"  Total unique genes: {n_total}\n"
            f"  Recognized:         {n_recognized} ({pct:.1f}%)\n"
            f"  Aliases resolved:   {len(aliases_resolved)}\n"
            f"  Unrecognized:       {len(unrecognized_symbols)}",
            err=True,
        )
        if unrecognized_symbols:
            listed = ", ".join(unrecognized_symbols[:10])
            click.echo(
                f"\n  Warning: {len(unrecognized_symbols)} unrecognized "
                f"symbols will be included as-is ({listed}). They will "
                "not match any gene in downstream analyses.",
                err=True,
            )

    # 3. Filter by --min-genes.
    if min_genes > 0:
        before = len(gene_sets)
        removed = [name for name, genes in gene_sets.items() if len(genes) < min_genes]
        gene_sets = {
            name: genes for name, genes in gene_sets.items() if len(genes) >= min_genes
        }
        if removed:
            click.echo(
                f"\nRemoved {len(removed)} gene sets with < {min_genes} "
                f"genes: {', '.join(removed)}",
                err=True,
            )

    # 4. Background genes.
    background_genes = None
    if background:
        background_genes = load_gene_set(path=background)
        click.echo(
            f"Background: {len(background_genes)} genes from {background}",
            err=True,
        )

    # 5. Build GeneSetCollection.
    collection = load_gene_sets_from_dict(
        gene_sets_dict=gene_sets,
        background_genes=background_genes,
        source="prepare-geneset",
    )

    # Attach metadata.
    input_basename = os.path.basename(input_path)
    collection.source_description = f"Custom gene sets from {input_basename}"
    collection.metadata = {
        "created_by": "hvantk prepare-geneset",
        "input_file": input_basename,
        "hgnc_validated": hgnc_validated,
    }
    if hgnc_validated:
        collection.metadata["hgnc_path"] = hgnc
        if aliases_resolved:
            collection.metadata["aliases_resolved"] = aliases_resolved
        if unrecognized_symbols:
            collection.metadata["unrecognized_symbols"] = unrecognized_symbols

    # 6. Save outputs.
    output_path.parent.mkdir(parents=True, exist_ok=True)
    collection.save(output_path)
    click.echo(f"\nSaved to {output}", err=True)

    if export_gmt:
        gmt_path = Path(export_gmt)
        gmt_path.parent.mkdir(parents=True, exist_ok=True)
        collection.save_gmt(gmt_path)
        click.echo(f"Exported GMT to {export_gmt}", err=True)
