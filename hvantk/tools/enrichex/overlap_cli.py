"""
CLI command for overlap enrichment testing using Fisher's exact test.

This module provides the `overlap` command for testing whether a query
gene list is enriched in gene sets.
"""

import logging

import click

logger = logging.getLogger(__name__)


def register_overlap_commands(group):
    """Register overlap commands to the enrichex CLI group.

    Parameters
    ----------
    group : click.Group
        The enrichex CLI group
    """
    group.add_command(overlap_test)


@click.command(name="overlap")
@click.option(
    "-g",
    "--gene-list",
    type=click.Path(exists=True),
    required=True,
    help="Path to query gene list (one gene per line, plain text)",
)
@click.option(
    "-s",
    "--gene-sets",
    type=click.Path(exists=True),
    required=True,
    help="Path to gene sets file (JSON format)",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    required=True,
    help="Output path for results",
)
@click.option(
    "--correction",
    type=click.Choice(["bonferroni", "benjamini-hochberg", "none"]),
    default="benjamini-hochberg",
    show_default=True,
    help="Multiple testing correction method",
)
@click.option(
    "--alpha",
    type=float,
    default=0.05,
    show_default=True,
    help="Significance threshold for reporting",
)
@click.option(
    "--output-format",
    type=click.Choice(["tsv", "json", "both"]),
    default="tsv",
    show_default=True,
    help="Output format",
)
@click.option(
    "--generate-report",
    is_flag=True,
    help="Generate plots and HTML report in the output directory",
)
@click.pass_context
def overlap_test(
    ctx, gene_list, gene_sets, output, correction, alpha, output_format, generate_report
):
    """Test gene list enrichment in gene sets using Fisher's exact test.

    This command tests whether a query gene list is significantly enriched
    in any of the provided gene sets using Fisher's exact test. Results
    include p-values, odds ratios, confidence intervals, and overlap genes.

    Gene sets should be prepared in advance using external tools or
    programmatically. The gene sets file should be in JSON format as
    produced by GeneSetCollection.save().

    \b
    Examples:
        # Basic usage
        hvantk enrichex overlap \\
            -g my_genes.txt \\
            -s cell_type_gene_sets.json \\
            -o enrichment_results.tsv

        # With Bonferroni correction
        hvantk enrichex overlap \\
            -g gwas_hits.txt \\
            -s gene_sets.json \\
            -o results.tsv \\
            --correction bonferroni

        # Output both TSV and JSON
        hvantk enrichex overlap \\
            -g genes.txt \\
            -s gene_sets.json \\
            -o results \\
            --output-format both

    \b
    Input file formats:
        Gene list: Plain text file with one gene per line
            Example:
                BRCA1
                TP53
                EGFR

        Gene sets: JSON file with GeneSetCollection format
            Created by: GeneSetCollection.save(path)
    """
    from hvantk.algorithms.enrichex.overlap import compute_overlap_enrichment_pandas
    from hvantk.core.utils.gene_sets import GeneSetCollection, load_gene_set

    # Load query genes
    click.echo(f"Loading query genes from: {gene_list}")
    query_genes = list(load_gene_set(path=gene_list))
    click.echo(f"  {len(query_genes)} genes loaded")

    # Load gene sets
    click.echo(f"Loading gene sets from: {gene_sets}")
    gene_set_collection = GeneSetCollection.load(gene_sets)
    click.echo(f"  {len(gene_set_collection)} gene sets loaded")
    click.echo(
        f"  Background universe: {len(gene_set_collection.background_genes)} genes"
    )

    # Run enrichment
    click.echo("\nRunning overlap enrichment analysis...")
    click.echo(f"  Correction method: {correction}")
    click.echo(f"  Significance threshold: {alpha}")

    df = compute_overlap_enrichment_pandas(query_genes, gene_set_collection, correction)

    if df.empty:
        click.echo("\nNo results generated (check gene ID matching)")
        return

    # Output results
    if output_format in ("tsv", "both"):
        tsv_path = output if output.endswith(".tsv") else f"{output}.tsv"
        df.to_csv(tsv_path, sep="\t", index=False)
        click.echo(f"\nTSV results written to: {tsv_path}")

    if output_format in ("json", "both"):
        json_path = output if output.endswith(".json") else f"{output}.json"
        df.to_json(json_path, orient="records", indent=2)
        click.echo(f"JSON results written to: {json_path}")

    # Generate report if requested
    if generate_report:
        from pathlib import Path

        from hvantk.algorithms.enrichex.report import generate_report

        click.echo("\nGenerating report...")

        # Determine output directory and file paths
        output_path = Path(output)
        if output_path.suffix in [".tsv", ".json"]:
            output_dir = output_path.parent
            results_path = output_path
        else:
            output_dir = output_path
            output_dir.mkdir(parents=True, exist_ok=True)
            results_path = output_dir / "overlap_results.tsv"
            if not results_path.exists():
                df.to_csv(results_path, sep="\t", index=False)

        report_path = output_dir / "enrichex_overlap_report.html"

        generate_report(
            output_path=str(report_path),
            overlap_results=str(results_path),
            gene_sets_path=gene_sets,
            title="EnrichEx Overlap Enrichment Report",
            description=f"Fisher's exact test for gene list enrichment in gene sets. Correction: {correction}, alpha: {alpha}",
            top_n=20,
            include_methods=True,
            include_gene_lists=True,
            embed_static_plots=True,
        )

        click.echo(f"✓ Report generated: {report_path}")
        click.echo(f"  Plots saved in: {output_dir}")

    # Summary statistics
    n_significant = (df["p_adjusted"] < alpha).sum()
    click.echo(f"\n{'='*60}")
    click.echo("ENRICHMENT SUMMARY")
    click.echo(f"{'='*60}")
    click.echo(f"Total gene sets tested: {len(df)}")
    click.echo(f"Significant gene sets (p_adj < {alpha}): {n_significant}")

    if n_significant > 0:
        click.echo("\nTop enriched gene sets:")
        click.echo("-" * 60)

        # Format top results table
        top = df.head(5)[
            ["gene_set_name", "n_overlap", "odds_ratio", "p_value", "p_adjusted"]
        ]
        # Format numbers for display
        top_display = top.copy()
        top_display["odds_ratio"] = top_display["odds_ratio"].apply(
            lambda x: f"{x:.2f}"
        )
        top_display["p_value"] = top_display["p_value"].apply(lambda x: f"{x:.2e}")
        top_display["p_adjusted"] = top_display["p_adjusted"].apply(
            lambda x: f"{x:.2e}"
        )

        # Print as table
        header = "  ".join(
            [
                "Gene Set".ljust(20),
                "Overlap".ljust(8),
                "OR".ljust(8),
                "P-value".ljust(10),
                "P-adj".ljust(10),
            ]
        )
        click.echo(header)
        click.echo("-" * 60)

        for _, row in top_display.iterrows():
            line = "  ".join(
                [
                    str(row["gene_set_name"])[:20].ljust(20),
                    str(row["n_overlap"]).ljust(8),
                    str(row["odds_ratio"]).ljust(8),
                    str(row["p_value"]).ljust(10),
                    str(row["p_adjusted"]).ljust(10),
                ]
            )
            click.echo(line)

    click.echo(f"{'='*60}\n")
