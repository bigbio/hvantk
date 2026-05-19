"""
Unified CLI for gene set extraction and preparation.

Subcommands:
  hvantk genesets clingen  -- Extract from ClinGen Hail Table (by GCEP/disease/keyword)
  hvantk genesets gencc    -- Extract from GenCC Hail Table (by submitter/disease/keyword)
  hvantk genesets prepare  -- Convert a plain TSV into a GeneSetCollection

The output is a GeneSetCollection file (JSON or GMT) consumable by
PSROC (``--gene-sets``), EnrichEx, or any compatible tool.
"""

import logging
import os
import sys
from pathlib import Path

import click

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------
# Shared output helpers
# ------------------------------------------------------------------


def _save_and_report(
    gene_sets,
    output_path,
    group_by,
    min_classification,
    min_genes,
    source,
    ctx,
):
    """Remove empty sets, report statistics, and save the collection."""
    from hvantk.core.utils.gene_sets import load_gene_sets_from_dict

    # Remove empty groups
    empty = [k for k, v in gene_sets.items() if not v]
    for k in empty:
        del gene_sets[k]
    if empty:
        click.echo(
            f"Warning: removed {len(empty)} empty group(s): "
            f"{', '.join(sorted(empty))}",
            err=True,
        )

    if not gene_sets:
        click.echo(
            "Error: No gene sets produced. " "Check data and filter settings.",
            err=True,
        )
        ctx.exit(1)

    # Report
    total_genes = len(set().union(*gene_sets.values()))
    sizes = sorted(len(g) for g in gene_sets.values())
    click.echo(
        f"Built {len(gene_sets)} gene sets "
        f"(group-by={group_by}, min-classification={min_classification}"
        + (f", min-genes={min_genes}" if min_genes > 0 else "")
        + ")"
    )
    click.echo(
        f"  Total unique genes: {total_genes}  |  "
        f"genes/group: min={sizes[0]}, "
        f"median={sizes[len(sizes) // 2]}, max={sizes[-1]}"
    )
    for name in sorted(gene_sets, key=lambda k: -len(gene_sets[k])):
        click.echo(f"  {name}: {len(gene_sets[name])} genes")

    # Save
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    collection = load_gene_sets_from_dict(
        {k: sorted(v) for k, v in gene_sets.items()},
        source=source,
    )
    if output_path.suffix.lower() == ".gmt":
        collection.save_gmt(str(output_path))
    else:
        collection.save(str(output_path))

    click.echo(f"\nSaved to: {output_path}")


def _check_overwrite(output_path, overwrite, ctx):
    """Exit with error if output exists and overwrite is False."""
    output_path = Path(output_path)
    if output_path.exists() and not overwrite:
        click.echo(
            f"Error: Output file already exists: {output_path}\n"
            "Use --overwrite to replace it.",
            err=True,
        )
        ctx.exit(1)


# ------------------------------------------------------------------
# Group definition
# ------------------------------------------------------------------


@click.group(name="genesets")
def genesets_group():
    """Extract or prepare gene set collections.

    \b
    Subcommands:
      clingen   Extract gene sets from ClinGen Hail Table
      gencc     Extract gene sets from GenCC Hail Table
      prepare   Convert a plain TSV into a GeneSetCollection
    """
    pass


# ------------------------------------------------------------------
# clingen subcommand
# ------------------------------------------------------------------


@genesets_group.command("clingen")
@click.option(
    "--ht",
    "clingen_ht",
    type=click.Path(exists=True),
    required=True,
    help="Path to ClinGen Hail Table (.ht) built with "
    "'hvantk mktable clingen-gene-disease'.",
)
@click.option(
    "--group-by",
    type=click.Choice(["gcep", "disease", "keyword"]),
    default="gcep",
    help="Grouping strategy [default: gcep]",
)
@click.option(
    "--min-classification",
    type=click.Choice(
        ["Definitive", "Strong", "Moderate", "Limited", "Disputed", "Refuted"]
    ),
    default="Moderate",
    help="Minimum classification level [default: Moderate]",
)
@click.option(
    "--min-genes",
    type=click.IntRange(min=0),
    default=0,
    help="Exclude groups with fewer than this many genes [default: 0]",
)
@click.option(
    "--categories-json",
    type=click.Path(exists=True),
    default=None,
    help="JSON file with keyword categories for --group-by keyword.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    required=True,
    help="Output path (.json or .gmt)",
)
@click.option("--overwrite", is_flag=True, help="Overwrite existing output file")
@click.pass_context
def genesets_clingen(
    ctx,
    clingen_ht,
    group_by,
    min_classification,
    min_genes,
    categories_json,
    output,
    overwrite,
):
    """Extract gene set collections from ClinGen gene-disease data.

    \b
    Grouping Strategies:
      gcep      Gene Curation Expert Panel (broader, recommended for PSROC)
      disease   Individual disease label (fine-grained)
      keyword   User-defined keyword categories (requires --categories-json)

    \b
    Examples:
      hvantk genesets clingen \\
          --ht /data/tables/clingen.ht \\
          --group-by gcep --min-genes 20 \\
          -o /data/gene_sets/clingen_gcep.json
    """
    import json

    _check_overwrite(output, overwrite, ctx)

    if group_by == "keyword" and not categories_json:
        click.echo(
            "Error: --categories-json is required when --group-by keyword",
            err=True,
        )
        ctx.exit(1)

    try:
        from hvantk.core.streamers.clingen import ClinGenStreamer

        streamer = ClinGenStreamer(clingen_ht)

        if group_by == "gcep":
            gene_sets = streamer.get_geneset_per_gcep(
                min_classification=min_classification,
                min_genes=min_genes,
            )
        elif group_by == "disease":
            gene_sets = streamer.get_geneset_per_disease(
                min_classification=min_classification,
            )
            if min_genes > 0:
                gene_sets = {k: v for k, v in gene_sets.items() if len(v) >= min_genes}
        else:  # keyword
            with open(categories_json, encoding="utf-8") as f:
                categories = json.load(f)
            gene_sets = streamer.aggregate_by_disease_category(
                categories=categories,
                min_classification=min_classification,
            )
            if min_genes > 0:
                gene_sets = {k: v for k, v in gene_sets.items() if len(v) >= min_genes}

        _save_and_report(
            gene_sets,
            output,
            group_by,
            min_classification,
            min_genes,
            f"clingen_{group_by}",
            ctx,
        )

    except click.exceptions.Exit:
        raise
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1) from e
    except Exception as e:
        logger.exception(f"Failed: {e}")
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1) from e


# ------------------------------------------------------------------
# gencc subcommand
# ------------------------------------------------------------------


@genesets_group.command("gencc")
@click.option(
    "--ht",
    "gencc_ht",
    type=click.Path(exists=True),
    required=True,
    help="Path to GenCC Hail Table (.ht) built with "
    "'hvantk mktable gencc-submissions'.",
)
@click.option(
    "--group-by",
    type=click.Choice(["submitter", "disease", "keyword"]),
    default="submitter",
    help="Grouping strategy [default: submitter]",
)
@click.option(
    "--min-classification",
    type=click.Choice(
        [
            "Definitive",
            "Strong",
            "Moderate",
            "Supportive",
            "Limited",
            "Disputed Evidence",
            "Refuted Evidence",
        ]
    ),
    default="Moderate",
    help="Minimum classification level [default: Moderate]",
)
@click.option(
    "--min-genes",
    type=click.IntRange(min=0),
    default=0,
    help="Exclude groups with fewer than this many genes [default: 0]",
)
@click.option(
    "--categories-json",
    type=click.Path(exists=True),
    default=None,
    help="JSON file with keyword categories for --group-by keyword.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    required=True,
    help="Output path (.json or .gmt)",
)
@click.option("--overwrite", is_flag=True, help="Overwrite existing output file")
@click.pass_context
def genesets_gencc(
    ctx,
    gencc_ht,
    group_by,
    min_classification,
    min_genes,
    categories_json,
    output,
    overwrite,
):
    """Extract gene set collections from GenCC submissions data.

    \b
    Grouping Strategies:
      submitter  Submitting organization (default)
      disease    Individual disease label
      keyword    User-defined keyword categories (requires --categories-json)

    \b
    Examples:
      hvantk genesets gencc \\
          --ht /data/tables/gencc.ht \\
          --group-by submitter --min-genes 20 \\
          -o /data/gene_sets/gencc_submitter.json
    """
    import json

    _check_overwrite(output, overwrite, ctx)

    if group_by == "keyword" and not categories_json:
        click.echo(
            "Error: --categories-json is required when --group-by keyword",
            err=True,
        )
        ctx.exit(1)

    try:
        from hvantk.core.streamers.gencc import GenCCStreamer

        streamer = GenCCStreamer(gencc_ht)

        if group_by == "submitter":
            gene_sets = streamer.get_geneset_per_submitter(
                min_classification=min_classification,
                min_genes=min_genes,
            )
        elif group_by == "disease":
            gene_sets = streamer.get_geneset_per_disease(
                min_classification=min_classification,
            )
            if min_genes > 0:
                gene_sets = {k: v for k, v in gene_sets.items() if len(v) >= min_genes}
        else:  # keyword
            with open(categories_json, encoding="utf-8") as f:
                categories = json.load(f)
            gene_sets = streamer.aggregate_by_disease_category(
                categories=categories,
                min_classification=min_classification,
            )
            if min_genes > 0:
                gene_sets = {k: v for k, v in gene_sets.items() if len(v) >= min_genes}

        _save_and_report(
            gene_sets,
            output,
            group_by,
            min_classification,
            min_genes,
            f"gencc_{group_by}",
            ctx,
        )

    except click.exceptions.Exit:
        raise
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1) from e
    except Exception as e:
        logger.exception(f"Failed: {e}")
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1) from e


# ------------------------------------------------------------------
# cosmic subcommand
# ------------------------------------------------------------------


@genesets_group.command("cosmic")
@click.option(
    "--ht",
    "cosmic_ht",
    type=click.Path(exists=True),
    required=True,
    help="Path to COSMIC CGC Hail Table (.ht) built with "
    "'hvantk mktable cosmic-cgc'.",
)
@click.option(
    "--group-by",
    type=click.Choice(["tumour-type", "role", "tissue", "keyword"]),
    default="tumour-type",
    help="Grouping strategy [default: tumour-type]",
)
@click.option(
    "--min-classification",
    type=click.Choice(["Tier 1", "Tier 2"]),
    default="Tier 1",
    help="Minimum classification level [default: Tier 1]",
)
@click.option(
    "--mutation-context",
    type=click.Choice(["somatic", "germline", "both"]),
    default="both",
    help="Filter genes by somatic/germline mutation context [default: both]",
)
@click.option(
    "--min-genes",
    type=click.IntRange(min=0),
    default=0,
    help="Exclude groups with fewer than this many genes [default: 0]",
)
@click.option(
    "--categories-json",
    type=click.Path(exists=True),
    default=None,
    help="JSON file with keyword categories for --group-by keyword.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    required=True,
    help="Output path (.json or .gmt)",
)
@click.option("--overwrite", is_flag=True, help="Overwrite existing output file")
@click.pass_context
def genesets_cosmic(
    ctx,
    cosmic_ht,
    group_by,
    min_classification,
    mutation_context,
    min_genes,
    categories_json,
    output,
    overwrite,
):
    """Extract gene set collections from COSMIC Cancer Gene Census data.

    \b
    Grouping Strategies:
      tumour-type  Tumour type (explodes multi-value field, default)
      role         Role in cancer (oncogene, TSG, fusion)
      tissue       Tissue type (E, L, M, O)
      keyword      User-defined keyword categories (requires --categories-json)

    \b
    Examples:
      hvantk genesets cosmic \\
          --ht /data/tables/cosmic_cgc.ht \\
          --group-by tumour-type \\
          --mutation-context somatic \\
          --min-classification "Tier 1" \\
          -o /data/gene_sets/cosmic_somatic_by_tumour.json
    """
    import json

    _check_overwrite(output, overwrite, ctx)

    if group_by == "keyword" and not categories_json:
        click.echo(
            "Error: --categories-json is required when --group-by keyword",
            err=True,
        )
        ctx.exit(1)

    try:
        from hvantk.core.streamers.cosmic_cgc import CosmicCGCStreamer

        streamer = CosmicCGCStreamer(cosmic_ht)

        if group_by == "tumour-type":
            gene_sets = streamer.get_geneset_per_tumour_type(
                min_classification=min_classification,
                min_genes=min_genes,
                mutation_context=mutation_context,
            )
        elif group_by == "role":
            gene_sets = streamer.get_geneset_per_role(
                min_classification=min_classification,
                min_genes=min_genes,
                mutation_context=mutation_context,
            )
        elif group_by == "tissue":
            gene_sets = streamer.get_geneset_per_tissue(
                min_classification=min_classification,
                min_genes=min_genes,
                mutation_context=mutation_context,
            )
        else:  # keyword
            with open(categories_json, encoding="utf-8") as f:
                categories = json.load(f)
            gene_sets = streamer.aggregate_by_disease_category(
                categories=categories,
                min_classification=min_classification,
            )
            if min_genes > 0:
                gene_sets = {k: v for k, v in gene_sets.items() if len(v) >= min_genes}

        _save_and_report(
            gene_sets,
            output,
            group_by,
            min_classification,
            min_genes,
            f"cosmic_{group_by}",
            ctx,
        )

    except click.exceptions.Exit:
        raise
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1) from e
    except Exception as e:
        logger.exception(f"Failed: {e}")
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1) from e


# ------------------------------------------------------------------
# prepare subcommand (migrated from prepare_geneset_cli.py)
# ------------------------------------------------------------------


@genesets_group.command("prepare")
@click.option(
    "-i",
    "--input",
    "input_path",
    type=click.Path(exists=True),
    required=True,
    help="Input TSV file (headerless, two columns: gene_set_name<TAB>gene_symbol).",
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
    help="Path to HGNC data (.ht or .tsv) for symbol validation and alias resolution.",
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
def genesets_prepare(
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
      hvantk genesets prepare -i panels.tsv -o panels.json --hgnc /data/hgnc.ht
    """
    from hvantk.core.utils.gene_sets import load_gene_set, load_gene_sets_from_dict
    from hvantk.core.utils.geneset_io import parse_geneset_tsv, validate_with_hgnc

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
        "created_by": "hvantk genesets prepare",
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
