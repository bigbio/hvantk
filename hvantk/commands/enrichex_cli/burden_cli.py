"""
CLI command for burden testing using Hail-native regression.

This module provides the `burden` command for testing whether cases have
excess rare variants in gene set genes using Hail's native logistic or
linear regression functions.
"""

import logging

import click

logger = logging.getLogger(__name__)


def register_burden_commands(group):
    """Register burden commands to the enrichex CLI group.

    Parameters
    ----------
    group : click.Group
        The enrichex CLI group
    """
    group.add_command(burden_test)


@click.command(name="burden")
@click.option(
    "-m",
    "--cohort-mt",
    type=click.Path(exists=True),
    required=True,
    help="Path to cohort MatrixTable (with variant annotations)",
)
@click.option(
    "-p",
    "--phenotypes",
    type=click.Path(exists=True),
    required=True,
    help="Path to phenotype file (Hail Table or TSV)",
)
@click.option(
    "-s",
    "--gene-sets",
    type=click.Path(exists=True),
    required=True,
    help="Path to gene sets file (JSON)",
)
@click.option(
    "-o", "--output", type=click.Path(), required=True, help="Output path for results"
)
@click.option(
    "--phenotype-field",
    type=str,
    default="is_case",
    show_default=True,
    help="Field name for phenotype",
)
@click.option(
    "--phenotype-type",
    type=click.Choice(["binary", "continuous"]),
    default="binary",
    show_default=True,
    help="Type of phenotype",
)
@click.option(
    "--covariates",
    type=str,
    default=None,
    help='Comma-separated covariate field names (e.g., "PC1,PC2,sex")',
)
@click.option(
    "--genotype-aggregation",
    type=click.Choice(["hets", "homs", "chets", "homs_chets"]),
    default="hets",
    show_default=True,
    help="Genotype aggregation method",
)
@click.option(
    "--max-af",
    type=float,
    default=0.01,
    show_default=True,
    help="Maximum allele frequency for qualifying variants",
)
@click.option(
    "--min-cadd",
    type=float,
    default=20,
    show_default=True,
    help="Minimum CADD score for qualifying variants (0 to disable)",
)
@click.option(
    "--consequences",
    type=str,
    default=None,
    help="Comma-separated qualifying consequences",
)
@click.option(
    "--gene-field",
    type=str,
    default="SYMBOL",
    show_default=True,
    help="Row field containing gene symbol",
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
    help="Significance threshold",
)
@click.option(
    "--sample-id-field",
    type=str,
    default="sample_id",
    show_default=True,
    help="Column name for sample ID in phenotype file (TSV only)",
)
@click.option("--dry-run", is_flag=True, help="Show execution plan without running")
@click.option(
    "--generate-report",
    is_flag=True,
    help="Generate plots and HTML report in the output directory",
)
@click.pass_context
def burden_test(
    ctx,
    cohort_mt,
    phenotypes,
    gene_sets,
    output,
    phenotype_field,
    phenotype_type,
    covariates,
    genotype_aggregation,
    max_af,
    min_cadd,
    consequences,
    gene_field,
    correction,
    alpha,
    sample_id_field,
    dry_run,
    generate_report,
):
    """Run case-control burden analysis using Hail-native regression.

    This command tests whether cases have excess rare damaging variants in
    gene set genes using Hail's native logistic_regression_rows() or
    linear_regression_rows() functions. The implementation follows the
    pattern from scripts/logreg_burden_test.py.

    Gene sets should be prepared in advance using external tools or
    programmatically. The gene sets file should be in JSON format as
    produced by GeneSetCollection.save().

    \b
    Workflow:
    1. Filter variants by AF, CADD, and consequence
    2. Aggregate variants to genes per sample (hets/homs/compound hets)
    3. Aggregate genes to gene sets per sample (burden score)
    4. Run logistic or linear regression with covariates
    5. Apply multiple testing correction

    \b
    Examples:
        # Basic logistic regression (binary phenotype)
        hvantk enrichex burden \\
            -m cohort.mt \\
            -p phenotypes.ht \\
            -s gene_sets.json \\
            -o burden_results.tsv

        # With covariates
        hvantk enrichex burden \\
            -m cohort.mt \\
            -p phenotypes.ht \\
            -s gene_sets.json \\
            --covariates PC1,PC2,PC3,PC4,PC5,sex \\
            -o burden_results.tsv

        # Linear regression (continuous phenotype)
        hvantk enrichex burden \\
            -m cohort.mt \\
            -p phenotypes.ht \\
            -s gene_sets.json \\
            --phenotype-type continuous \\
            --phenotype-field cognitive_score \\
            -o burden_results.tsv

        # Recessive model (homs + compound hets)
        hvantk enrichex burden \\
            -m cohort.mt \\
            -p phenotypes.ht \\
            -s gene_sets.json \\
            --genotype-aggregation homs_chets \\
            --max-af 0.001 \\
            -o burden_results.tsv

        # From TSV phenotype file
        hvantk enrichex burden \\
            -m cohort.mt \\
            -p phenotypes.tsv \\
            -s gene_sets.json \\
            --sample-id-field sample_id \\
            -o burden_results.tsv

    \b
    Input Requirements:
        cohort.mt: Hail MatrixTable with:
        - Row fields: gene annotations (SYMBOL), variant annotations (gnomad_af, cadd_phred)
        - Entry fields: genotypes (GT), quality (GQ, DP)

        phenotypes: Hail Table or TSV with:
        - Sample IDs (must match MT column keys)
        - Phenotype field (binary or continuous)
        - Optional covariate fields

        gene_sets.json: JSON file created by GeneSetCollection.save()
    """
    from hvantk.core.hail_context import init_hail

    # Initialize Hail
    init_hail()

    # Show plan if dry-run
    if dry_run:
        click.echo("\n" + "=" * 60)
        click.echo("ENRICHEX BURDEN ANALYSIS PLAN (Hail-Native)")
        click.echo("=" * 60)
        click.echo(f"\nInputs:")
        click.echo(f"  Cohort MT: {cohort_mt}")
        click.echo(f"  Phenotypes: {phenotypes}")
        click.echo(f"    Field: {phenotype_field}")
        click.echo(f"    Type: {phenotype_type}")
        click.echo(f"  Gene sets: {gene_sets}")
        click.echo(f"\nVariant Filters:")
        click.echo(f"  Max AF: {max_af}")
        click.echo(f"  Min CADD: {min_cadd if min_cadd > 0 else 'disabled'}")
        if consequences:
            click.echo(f"  Consequences: {consequences}")
        click.echo(f"\nAnalysis:")
        click.echo(f"  Gene field: {gene_field}")
        click.echo(f"  Genotype aggregation: {genotype_aggregation}")
        if covariates:
            click.echo(f"  Covariates: {covariates}")
        click.echo(f"  Correction: {correction}")
        click.echo(f"\nOutput: {output}")
        click.echo("=" * 60 + "\n")
        return

    import hail as hl

    from hvantk.enrichex.burden import run_burden_analysis
    from hvantk.enrichex.correction import apply_correction
    from hvantk.enrichex.gene_sets import GeneSetCollection

    # Load inputs
    click.echo(f"Loading cohort: {cohort_mt}")
    mt = hl.read_matrix_table(cohort_mt)
    click.echo(f"  Rows: {mt.count_rows()}, Cols: {mt.count_cols()}")

    click.echo(f"\nLoading phenotypes: {phenotypes}")
    if phenotypes.endswith(".ht"):
        phenotype_ht = hl.read_table(phenotypes)
    else:
        # Import TSV and convert to Hail Table
        phenotype_ht = hl.import_table(phenotypes, impute=True).key_by(sample_id_field)
    click.echo(f"  {phenotype_ht.count()} samples with phenotypes")

    click.echo(f"\nLoading gene sets: {gene_sets}")
    gene_set_collection = GeneSetCollection.load(gene_sets)
    click.echo(f"  {len(gene_set_collection)} gene sets loaded")

    # Convert to dict for burden function
    gene_sets_dict = {gs.name: list(gs.genes) for gs in gene_set_collection}

    # Parse options
    conseq_list = consequences.split(",") if consequences else None
    covar_list = covariates.split(",") if covariates else None

    # Disable CADD filter if set to 0
    cadd_threshold = min_cadd if min_cadd > 0 else None

    # Run Hail-native burden analysis
    click.echo("\n" + "=" * 60)
    click.echo("RUNNING BURDEN ANALYSIS")
    click.echo("=" * 60 + "\n")

    result_ht = run_burden_analysis(
        cohort_mt=mt,
        gene_sets=gene_sets_dict,
        phenotype_ht=phenotype_ht,
        phenotype_field=phenotype_field,
        covariate_fields=covar_list,
        phenotype_type=phenotype_type,
        gene_field=gene_field,
        max_af=max_af,
        min_cadd=cadd_threshold,
        consequences=conseq_list,
        genotype_aggregation=genotype_aggregation,
    )

    # Convert to pandas for easier manipulation
    click.echo("\nConverting results to pandas...")
    result_df = result_ht.to_pandas()

    # Apply multiple testing correction
    click.echo(f"Applying {correction} correction...")
    p_values = result_df["p_value"].tolist()
    p_adjusted = apply_correction(p_values, method=correction)
    result_df["p_adjusted"] = p_adjusted
    result_df["significant"] = result_df["p_adjusted"] < alpha

    # Sort by p-value
    result_df = result_df.sort_values("p_value")

    # Output (resolve output path once to allow directory targets)
    from pathlib import Path

    output_path = Path(output)
    if output_path.suffix in [".tsv", ".json"]:
        output_dir = output_path.parent
        results_path = output_path
    else:
        output_dir = output_path
        output_dir.mkdir(parents=True, exist_ok=True)
        results_path = output_dir / "burden_results.tsv"

    # Output
    result_df.to_csv(results_path, sep="\t", index=False)
    click.echo(f"\n✓ Results written to: {results_path}")

    # Summary statistics
    n_significant = result_df["significant"].sum()
    click.echo("\n" + "=" * 60)
    click.echo("BURDEN ANALYSIS SUMMARY")
    click.echo("=" * 60)
    click.echo(f"Total gene sets tested: {len(result_df)}")
    click.echo(f"Significant gene sets (p_adj < {alpha}): {n_significant}")

    if n_significant > 0:
        click.echo(f"\nTop significant gene sets:")
        click.echo("-" * 60)

        # Select columns based on phenotype type
        if phenotype_type == "binary":
            display_cols = [
                "gene_set_name",
                "beta",
                "odds_ratio",
                "p_value",
                "p_adjusted",
            ]
        else:
            display_cols = [
                "gene_set_name",
                "beta",
                "standard_error",
                "p_value",
                "p_adjusted",
            ]

        # Filter to available columns
        available_cols = [c for c in display_cols if c in result_df.columns]
        top = result_df[result_df["significant"]].head(5)[available_cols]

        # Format for display
        top_display = top.copy()
        if "beta" in top_display.columns:
            top_display["beta"] = top_display["beta"].apply(lambda x: f"{x:.4f}")
        if "odds_ratio" in top_display.columns:
            top_display["odds_ratio"] = top_display["odds_ratio"].apply(
                lambda x: f"{x:.2f}"
            )
        if "standard_error" in top_display.columns:
            top_display["standard_error"] = top_display["standard_error"].apply(
                lambda x: f"{x:.4f}"
            )
        top_display["p_value"] = top_display["p_value"].apply(lambda x: f"{x:.2e}")
        top_display["p_adjusted"] = top_display["p_adjusted"].apply(
            lambda x: f"{x:.2e}"
        )

        # Print table
        if phenotype_type == "binary":
            header = "  ".join(
                [
                    "Gene Set".ljust(20),
                    "Beta".ljust(8),
                    "OR".ljust(8),
                    "P-value".ljust(10),
                    "P-adj".ljust(10),
                ]
            )
        else:
            header = "  ".join(
                [
                    "Gene Set".ljust(20),
                    "Beta".ljust(8),
                    "SE".ljust(8),
                    "P-value".ljust(10),
                    "P-adj".ljust(10),
                ]
            )

        click.echo(header)
        click.echo("-" * 60)

        for _, row in top_display.iterrows():
            if phenotype_type == "binary":
                line = "  ".join(
                    [
                        str(row["gene_set_name"])[:20].ljust(20),
                        str(row["beta"]).ljust(8),
                        str(row.get("odds_ratio", "N/A")).ljust(8),
                        str(row["p_value"]).ljust(10),
                        str(row["p_adjusted"]).ljust(10),
                    ]
                )
            else:
                line = "  ".join(
                    [
                        str(row["gene_set_name"])[:20].ljust(20),
                        str(row["beta"]).ljust(8),
                        str(row.get("standard_error", "N/A")).ljust(8),
                        str(row["p_value"]).ljust(10),
                        str(row["p_adjusted"]).ljust(10),
                    ]
                )
            click.echo(line)

    click.echo("=" * 60 + "\n")

    # Generate report if requested
    if generate_report:
        from hvantk.enrichex.report import generate_report as gen_report

        click.echo("\nGenerating report...")

        report_path = output_dir / "enrichex_burden_report.html"

        gen_report(
            output_path=report_path,
            burden_results=results_path,
            gene_sets_path=gene_sets,
            phenotype_type=phenotype_type,
            top_n=25,
            embed_static_plots=True,
        )
        click.echo(f"\n✓ Report written to: {report_path}")
