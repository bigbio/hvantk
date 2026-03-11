"""
ClinGen gene set extraction CLI.

Extracts gene set collections from ClinGen gene-disease validity data,
grouped by GCEP, disease label, or user-defined keyword categories.

The output is a GeneSetCollection file (JSON or GMT) consumable by
PSROC (``--gene-sets``), EnrichEx, or any compatible tool.
"""

import logging
import click

logger = logging.getLogger(__name__)


@click.command(name="clingen-genesets")
@click.option(
    "--clingen-ht",
    type=click.Path(exists=True),
    required=True,
    help="Path to ClinGen Hail Table (.ht) built with "
    "'hvantk mktable clingen-gene-disease'.",
)
@click.option(
    "--group-by",
    type=click.Choice(["gcep", "disease", "keyword"]),
    default="gcep",
    help="Grouping strategy: 'gcep' (expert panel, recommended), "
    "'disease' (individual disease label), "
    "'keyword' (custom categories via --categories-json) [default: gcep]",
)
@click.option(
    "--min-classification",
    type=click.Choice(
        [
            "Definitive",
            "Strong",
            "Moderate",
            "Limited",
            "Disputed",
            "Refuted",
        ]
    ),
    default="Moderate",
    help="Minimum ClinGen classification level [default: Moderate]",
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
    help="JSON file with keyword categories for --group-by keyword. "
    'Format: {"name": ["term1", "term2", ...]}',
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    required=True,
    help="Output path (.json or .gmt)",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite existing output file",
)
@click.pass_context
def clingen_genesets_cmd(
    ctx,
    clingen_ht,
    group_by,
    min_classification,
    min_genes,
    categories_json,
    output,
    overwrite,
):
    """
    Extract gene set collections from ClinGen gene-disease data.

    \b
    Grouping Strategies:
      gcep      Gene Curation Expert Panel (broader, recommended for PSROC)
      disease   Individual disease label (fine-grained, most have 1-3 genes)
      keyword   User-defined keyword categories (requires --categories-json)

    \b
    Examples:

      # GCEP grouping with minimum 20 genes per panel
      hvantk clingen-genesets \\
          --clingen-ht /data/tables/clingen.ht \\
          --group-by gcep \\
          --min-classification Moderate \\
          --min-genes 20 \\
          -o /data/gene_sets/clingen_gcep.json

      # Then feed into PSROC:
      hvantk psroc \\
          --gene-sets /data/gene_sets/clingen_gcep.json \\
          --clinvar-ht /data/clinvar.ht \\
          --dbnsfp-ht /data/dbnsfp.ht \\
          --scores "CADD_phred,REVEL_score" \\
          --output-dir /results/psroc_by_gcep
    """
    import json
    from pathlib import Path

    output_path = Path(output)
    if output_path.exists() and not overwrite:
        click.echo(
            f"Error: Output file already exists: {output_path}\n"
            "Use --overwrite to replace it.",
            err=True,
        )
        ctx.exit(1)

    if group_by == "keyword" and not categories_json:
        click.echo(
            "Error: --categories-json is required when --group-by keyword",
            err=True,
        )
        ctx.exit(1)

    try:
        from hvantk.data.clingen_streamer import ClinGenStreamer
        from hvantk.utils.gene_sets import load_gene_sets_from_dict

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
            with open(categories_json) as f:
                categories = json.load(f)
            gene_sets = streamer.aggregate_by_disease_category(
                categories=categories,
                min_classification=min_classification,
            )
            if min_genes > 0:
                gene_sets = {k: v for k, v in gene_sets.items() if len(v) >= min_genes}

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
                "Error: No gene sets produced. "
                "Check ClinGen data and filter settings.",
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
            f"median={sizes[len(sizes)//2]}, max={sizes[-1]}"
        )
        for name in sorted(gene_sets, key=lambda k: -len(gene_sets[k])):
            click.echo(f"  {name}: {len(gene_sets[name])} genes")

        # Save
        output_path.parent.mkdir(parents=True, exist_ok=True)
        collection = load_gene_sets_from_dict(
            {k: sorted(v) for k, v in gene_sets.items()},
            source=f"clingen_{group_by}",
        )
        if output_path.suffix.lower() == ".gmt":
            collection.save_gmt(str(output_path))
        else:
            collection.save(str(output_path))

        click.echo(f"\nSaved to: {output_path}")

    except click.exceptions.Exit:
        raise
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1)
    except Exception as e:
        logger.exception(f"Failed: {e}")
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1)
