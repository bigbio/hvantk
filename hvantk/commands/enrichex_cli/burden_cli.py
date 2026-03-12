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
    group.add_command(burden_pipeline_cmd)


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
    type=click.Choice(["hets", "homs", "multi_het", "homs_multi_het", "chets", "homs_chets"]),
    default="hets",
    show_default=True,
    help="Genotype aggregation method. 'multi_het' counts genes with >= 2 het variants "
    "(approximates compound-het). 'chets'/'homs_chets' are deprecated aliases.",
)
@click.option(
    "--max-af",
    type=float,
    default=None,
    help="Maximum allele frequency (deprecated, use VariantFilter API instead)",
)
@click.option(
    "--min-score",
    type=float,
    default=None,
    help="Minimum prediction score threshold (e.g., CADD, REVEL). None to disable.",
)
@click.option(
    "--consequences",
    type=str,
    default=None,
    help="Comma-separated qualifying consequences for variant-class selection",
)
@click.option(
    "--gene-field",
    type=str,
    default="SYMBOL",
    show_default=True,
    help="Row field containing gene symbol",
)
@click.option(
    "--af-field",
    type=str,
    default="gnomad_af",
    show_default=True,
    help="Row field containing allele frequency",
)
@click.option(
    "--score-field",
    type=str,
    default="cadd_phred",
    show_default=True,
    help="Row field containing prediction score (e.g., cadd_phred, REVEL)",
)
@click.option(
    "--consequence-field",
    type=str,
    default="consequence",
    show_default=True,
    help="Row field containing variant consequence",
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
@click.option(
    "--variant-classes",
    type=str,
    default=None,
    help="Comma-separated variant class presets for stratified analysis "
    "(e.g., 'lof,missense_constrained,synonymous'). When specified, burden "
    "is run independently for each class and outputs per-class TSVs plus "
    "a combined summary. Available presets: lof, missense_constrained, synonymous.",
)
@click.option(
    "--permutation",
    is_flag=True,
    help="Run permutation-based competitive gene-set test via gene-label permutation "
    "in addition to the self-contained test.",
)
@click.option(
    "--n-permutations",
    type=int,
    default=10000,
    show_default=True,
    help="Number of permutations for permutation test.",
)
@click.option(
    "--permutation-seed",
    type=int,
    default=None,
    help="Random seed for permutation test reproducibility.",
)
@click.option(
    "--normalize-by-length",
    is_flag=True,
    help="Normalize per-gene burden by gene length before aggregation. "
    "Controls for longer genes accumulating more rare variants by chance.",
)
@click.option(
    "--gene-lengths",
    type=click.Path(exists=True),
    default=None,
    help="Path to TSV file with gene lengths (columns: gene, length_bp). "
    "Used with --normalize-by-length. If omitted, variant site count "
    "per gene is used as proxy.",
)
@click.option(
    "--min-carriers",
    type=int,
    default=0,
    show_default=True,
    help="Minimum number of samples carrying qualifying variants for a gene "
    "set to be included in regression. Gene sets with fewer carriers are "
    "skipped (regression would be uninformative).",
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
    min_score,
    consequences,
    gene_field,
    af_field,
    score_field,
    consequence_field,
    correction,
    alpha,
    sample_id_field,
    variant_classes,
    permutation,
    n_permutations,
    permutation_seed,
    normalize_by_length,
    gene_lengths,
    min_carriers,
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
    1. Optionally filter variants by consequence class (MT expected pre-filtered)
    2. Aggregate variants to genes per sample (hets/homs/multi-het)
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
        if max_af is not None:
            click.echo(f"  Max AF: {max_af} (field: {af_field})")
        if min_score is not None:
            click.echo(f"  Min score: {min_score if min_score > 0 else 'disabled'} (field: {score_field})")
        if consequences:
            click.echo(f"  Consequences: {consequences} (field: {consequence_field})")
        if max_af is None and min_score is None and not consequences:
            click.echo(f"  None (MT assumed pre-filtered)")
        click.echo(f"\nAnalysis:")
        click.echo(f"  Gene field: {gene_field}")
        click.echo(f"  Genotype aggregation: {genotype_aggregation}")
        if covariates:
            click.echo(f"  Covariates: {covariates}")
        if normalize_by_length:
            click.echo(f"  Normalize by length: Yes")
            if gene_lengths:
                click.echo(f"  Gene lengths file: {gene_lengths}")
            else:
                click.echo(f"  Gene lengths: using variant site count proxy")
        if min_carriers > 0:
            click.echo(f"  Min carriers: {min_carriers}")
        if variant_classes:
            click.echo(f"  Variant classes: {variant_classes}")
        if permutation:
            click.echo(f"  Permutation test: Yes ({n_permutations} permutations)")
        click.echo(f"  Correction: {correction}")
        click.echo(f"\nOutput: {output}")
        click.echo("=" * 60 + "\n")
        return

    import hail as hl
    import pandas as pd

    from hvantk.enrichex.burden import VariantFilter, run_burden_analysis
    from hvantk.enrichex.correction import apply_correction
    from hvantk.utils.gene_sets import GeneSetCollection

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

    # Load gene lengths if provided
    gene_lengths_dict = None
    if gene_lengths:
        click.echo(f"\nLoading gene lengths: {gene_lengths}")
        gene_lengths_df = pd.read_csv(gene_lengths, sep="\t")
        gene_lengths_dict = dict(
            zip(gene_lengths_df.iloc[:, 0], gene_lengths_df.iloc[:, 1])
        )
        click.echo(f"  {len(gene_lengths_dict)} gene lengths loaded")

    # Parse options
    conseq_list = consequences.split(",") if consequences else None
    covar_list = covariates.split(",") if covariates else None

    # Disable score filter if set to 0
    score_threshold = min_score if (min_score is not None and min_score > 0) else None

    # Resolve output path
    from pathlib import Path

    output_path = Path(output)
    if output_path.suffix in [".tsv", ".json"]:
        output_dir = output_path.parent
        results_path = output_path
    else:
        output_dir = output_path
        output_dir.mkdir(parents=True, exist_ok=True)
        results_path = output_dir / "burden_results.tsv"

    # --- Stratified analysis ---
    if variant_classes:
        from hvantk.enrichex.burden import (
            build_variant_classes_from_presets,
            run_stratified_burden_analysis,
        )

        class_names = [c.strip() for c in variant_classes.split(",")]

        # Build base filter from CLI params (for AF/field settings)
        has_filter_params = (
            max_af is not None or score_threshold is not None or conseq_list is not None
        )
        base_filter = None
        if has_filter_params:
            base_filter = VariantFilter(
                max_af=max_af if max_af is not None else 1.0,
                min_score=score_threshold,
                consequences=None,  # overridden per class
                pass_only=False,
                min_gq=0,
                min_dp=0,
                af_field=af_field,
                score_field=score_field,
                consequence_field=consequence_field,
            )
        else:
            # Use default field names from CLI options
            base_filter = VariantFilter(
                max_af=1.0,
                min_score=None,
                consequences=None,
                pass_only=False,
                min_gq=0,
                min_dp=0,
                af_field=af_field,
                score_field=score_field,
                consequence_field=consequence_field,
            )

        vc_dict = build_variant_classes_from_presets(class_names, base_filter=base_filter)

        click.echo("\n" + "=" * 60)
        click.echo("RUNNING STRATIFIED BURDEN ANALYSIS")
        click.echo("=" * 60)
        click.echo(f"Variant classes: {class_names}\n")

        stratified_results = run_stratified_burden_analysis(
            cohort_mt=mt,
            gene_sets=gene_sets_dict,
            phenotype_ht=phenotype_ht,
            variant_classes=vc_dict,
            phenotype_field=phenotype_field,
            covariate_fields=covar_list,
            phenotype_type=phenotype_type,
            gene_field=gene_field,
            genotype_aggregation=genotype_aggregation,
            normalize_by_length=normalize_by_length,
            gene_lengths=gene_lengths_dict,
            min_carriers=min_carriers,
        )

        if not stratified_results:
            click.echo(
                "\nNo variant classes produced testable results. "
                "Check that gene IDs match and samples overlap with phenotypes."
            )
            result_df = pd.DataFrame()
            result_df.to_csv(results_path, sep="\t", index=False)
            click.echo(f"\n  Empty results written to: {results_path}")
            return

        # Convert to pandas, apply correction, write per-class TSVs
        all_dfs = []
        for class_name, result_ht in stratified_results.items():
            df = result_ht.to_pandas()
            p_values = df["p_value"].tolist()
            p_adjusted = apply_correction(p_values, method=correction)
            df["p_adjusted"] = p_adjusted
            df["significant"] = df["p_adjusted"] < alpha
            df = df.sort_values("p_value")

            class_path = output_dir / f"burden_{class_name}.tsv"
            df.to_csv(class_path, sep="\t", index=False)
            click.echo(f"  {class_name}: {class_path}")
            all_dfs.append(df)

        # Combined output
        if all_dfs:
            combined_df = pd.concat(all_dfs, ignore_index=True)
            combined_df.to_csv(results_path, sep="\t", index=False)
            click.echo(f"\n  Combined: {results_path}")

        # Summary
        click.echo("\n" + "=" * 60)
        click.echo("STRATIFIED BURDEN ANALYSIS SUMMARY")
        click.echo("=" * 60)
        for class_name in stratified_results:
            class_df = next(df for df in all_dfs if df["variant_class"].iloc[0] == class_name)
            n_sig = class_df["significant"].sum()
            click.echo(
                f"  {class_name}: {len(class_df)} gene sets tested, "
                f"{n_sig} significant (p_adj < {alpha})"
            )
        click.echo("=" * 60 + "\n")

        result_df = combined_df if all_dfs else pd.DataFrame()

    # --- Single-class analysis ---
    else:
        # Build VariantFilter if any filter params provided
        has_filter_params = (
            max_af is not None or score_threshold is not None or conseq_list is not None
        )
        variant_filter_obj = None
        if has_filter_params:
            variant_filter_obj = VariantFilter(
                max_af=max_af if max_af is not None else 1.0,
                min_score=score_threshold,
                consequences=conseq_list,
                pass_only=False,
                min_gq=0,
                min_dp=0,
                af_field=af_field,
                score_field=score_field,
                consequence_field=consequence_field,
            )

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
            genotype_aggregation=genotype_aggregation,
            variant_filter=variant_filter_obj,
            normalize_by_length=normalize_by_length,
            gene_lengths=gene_lengths_dict,
            min_carriers=min_carriers,
        )

        if result_ht is None:
            click.echo(
                "\nNo testable gene sets — burden analysis produced no results. "
                "Check that gene IDs match the MatrixTable and that samples "
                "overlap with the phenotype table."
            )
            result_df = pd.DataFrame()
            result_df.to_csv(results_path, sep="\t", index=False)
            click.echo(f"\n  Empty results written to: {results_path}")
            return

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

        # Output
        result_df.to_csv(results_path, sep="\t", index=False)
        click.echo(f"\n  Results written to: {results_path}")

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

    # Run permutation test if requested
    if permutation:
        from hvantk.enrichex.burden import permutation_burden_test

        click.echo("\n" + "=" * 60)
        click.echo("RUNNING PERMUTATION BURDEN TEST")
        click.echo("=" * 60 + "\n")

        if variant_classes:
            # Stratified permutation: run for each variant class
            for class_name, vf in vc_dict.items():
                click.echo(f"\nPermutation test for variant class: {class_name}")
                comp_df = permutation_burden_test(
                    cohort_mt=mt,
                    gene_sets=gene_sets_dict,
                    phenotype_ht=phenotype_ht,
                    phenotype_field=phenotype_field,
                    covariate_fields=covar_list,
                    phenotype_type=phenotype_type,
                    gene_field=gene_field,
                    genotype_aggregation=genotype_aggregation,
                    variant_filter=vf,
                    n_permutations=n_permutations,
                    seed=permutation_seed,
                    normalize_by_length=normalize_by_length,
                    gene_lengths=gene_lengths_dict,
                )
                if not comp_df.empty:
                    comp_df["variant_class"] = class_name
                    comp_path = output_dir / f"permutation_{class_name}.tsv"
                    comp_df.to_csv(comp_path, sep="\t", index=False)
                    click.echo(f"    {comp_path}")
        else:
            # Single permutation test
            comp_df = permutation_burden_test(
                cohort_mt=mt,
                gene_sets=gene_sets_dict,
                phenotype_ht=phenotype_ht,
                phenotype_field=phenotype_field,
                covariate_fields=covar_list,
                phenotype_type=phenotype_type,
                gene_field=gene_field,
                genotype_aggregation=genotype_aggregation,
                variant_filter=variant_filter_obj if not variant_classes else None,
                n_permutations=n_permutations,
                seed=permutation_seed,
                normalize_by_length=normalize_by_length,
                gene_lengths=gene_lengths_dict,
            )
            if not comp_df.empty:
                comp_path = output_dir / "permutation_results.tsv"
                comp_df.to_csv(comp_path, sep="\t", index=False)
                click.echo(f"\n  Permutation results: {comp_path}")

    # Generate report if requested
    if generate_report:
        from hvantk.enrichex.report import generate_report as gen_report

        click.echo("\nGenerating report...")

        report_path = output_dir / "enrichex_burden_report.html"

        gen_report(
            output_path=report_path,
            burden_results=results_path,
            gene_sets_path=gene_sets,
            top_n=25,
            embed_static_plots=True,
        )
        click.echo(f"\n  Report written to: {report_path}")


# ---------------------------------------------------------------------------
# burden-pipeline subcommand
# ---------------------------------------------------------------------------


@click.command(name="burden-pipeline")
@click.option(
    "-m",
    "--cohort-mt",
    type=click.Path(exists=True),
    required=True,
    help="Path to cohort MatrixTable.",
)
@click.option(
    "-p",
    "--phenotypes",
    type=click.Path(exists=True),
    default=None,
    help="Path to phenotype file (Hail Table or TSV). "
    "If omitted, phenotype and covariates are extracted from MT column fields "
    "(use dot notation for nested structs, e.g., --phenotype-field phe.is_case).",
)
@click.option(
    "--gene-sets",
    "gene_sets_pairs",
    type=str,
    required=True,
    multiple=True,
    help="Gene set collection as name:path (repeatable). "
    "Example: --gene-sets heart:gene_sets/heart.json --gene-sets brain:gene_sets/brain.json",
)
@click.option(
    "--phenotype-field",
    type=str,
    default="is_case",
    show_default=True,
    help="Phenotype column name. Supports dot notation for nested MT column "
    "structs (e.g., 'phe.is_case').",
)
@click.option(
    "--phenotype-type",
    type=click.Choice(["binary", "continuous"]),
    default="binary",
    show_default=True,
    help="Type of phenotype.",
)
@click.option(
    "--covariates",
    type=str,
    default=None,
    help='Comma-separated covariate field names (e.g., "PC1,PC2,sex"). '
    "Supports dot notation for nested MT column structs.",
)
@click.option(
    "--variant-classes",
    type=str,
    default=None,
    help="Comma-separated variant class presets "
    "(e.g., 'lof,missense_constrained,synonymous'). "
    "If omitted, a single unstratified run is performed.",
)
@click.option(
    "--gene-field",
    type=str,
    default="SYMBOL",
    show_default=True,
    help="Row field containing gene symbol.",
)
@click.option(
    "--genotype-aggregation",
    type=click.Choice(
        ["hets", "homs", "multi_het", "homs_multi_het", "chets", "homs_chets"]
    ),
    default="hets",
    show_default=True,
    help="Genotype aggregation method.",
)
@click.option(
    "--correction",
    type=click.Choice(["bonferroni", "benjamini-hochberg", "none"]),
    default="benjamini-hochberg",
    show_default=True,
    help="Multiple testing correction method.",
)
@click.option(
    "--alpha",
    type=float,
    default=0.05,
    show_default=True,
    help="Significance threshold.",
)
@click.option(
    "--min-carriers",
    type=int,
    default=5,
    show_default=True,
    help="Minimum carriers per gene set for regression.",
)
@click.option(
    "--competitive",
    is_flag=True,
    help="Run permutation-based competitive testing.",
)
@click.option(
    "--n-permutations",
    type=int,
    default=10000,
    show_default=True,
    help="Number of permutations for competitive test.",
)
@click.option(
    "--permutation-seed",
    type=int,
    default=None,
    help="Random seed for permutation reproducibility.",
)
@click.option(
    "--normalize-by-length",
    is_flag=True,
    help="Normalize per-gene burden by gene length.",
)
@click.option(
    "--gene-lengths",
    type=click.Path(exists=True),
    default=None,
    help="TSV with gene lengths (columns: gene, length_bp).",
)
@click.option(
    "--af-field",
    type=str,
    default="gnomad_af",
    show_default=True,
    help="Row field containing allele frequency. Supports dot notation (e.g., 'gnomAD_AF').",
)
@click.option(
    "--score-field",
    type=str,
    default="cadd_phred",
    show_default=True,
    help="Row field containing prediction score (e.g., 'cadd_phred', 'REVEL', 'vep.CADD_PHRED'). "
    "Supports dot notation for nested structs.",
)
@click.option(
    "--consequence-field",
    type=str,
    default="consequence",
    show_default=True,
    help="Row field containing variant consequence. Supports dot notation (e.g., 'Consequence').",
)
@click.option(
    "--max-af",
    type=float,
    default=0.01,
    show_default=True,
    help="Maximum allele frequency threshold.",
)
@click.option(
    "--pass-only/--no-pass-only",
    default=False,
    show_default=True,
    help="Only include PASS variants (requires 'filters' field in MT).",
)
@click.option(
    "--sample-id-field",
    type=str,
    default="sample_id",
    show_default=True,
    help="Sample ID column name in TSV phenotype file.",
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    required=True,
    help="Output directory for all results.",
)
@click.option("--dry-run", is_flag=True, help="Show execution plan without running.")
@click.option(
    "--generate-report",
    is_flag=True,
    help="Generate HTML report after pipeline completes.",
)
def burden_pipeline_cmd(
    cohort_mt,
    phenotypes,
    gene_sets_pairs,
    phenotype_field,
    phenotype_type,
    covariates,
    variant_classes,
    gene_field,
    genotype_aggregation,
    correction,
    alpha,
    min_carriers,
    competitive,
    n_permutations,
    permutation_seed,
    normalize_by_length,
    gene_lengths,
    af_field,
    score_field,
    consequence_field,
    max_af,
    pass_only,
    sample_id_field,
    output_dir,
    dry_run,
    generate_report,
):
    """Run multi-tissue, multi-variant-class burden pipeline.

    This command orchestrates burden analysis across multiple gene set
    collections and variant classes in a single coordinated run.  It
    loads the cohort MatrixTable once and iterates over the cross-product
    of (variant_class x gene_set_collection).

    \b
    Examples:
        # Multi-tissue, multi-variant-class
        hvantk enrichex burden-pipeline \\
            -m cohort.mt \\
            -p phenotypes.ht \\
            --gene-sets heart:gene_sets/heart.json \\
            --gene-sets brain:gene_sets/brain.json \\
            --variant-classes lof,missense_constrained,synonymous \\
            --covariates PC1,PC2,PC3,sex \\
            --competitive \\
            --generate-report \\
            -o results/chd_celltype_burden/

        # Single tissue, no stratification
        hvantk enrichex burden-pipeline \\
            -m cohort.mt \\
            -p phenotypes.ht \\
            --gene-sets markers:gene_sets/markers.json \\
            -o results/burden/

        # Phenotype from MT column fields (no separate phenotype file)
        hvantk enrichex burden-pipeline \\
            -m cohort.mt \\
            --gene-sets heart:gene_sets/heart.json \\
            --phenotype-field phe.is_case \\
            --covariates phe.PC1,phe.PC2,phe.sex \\
            -o results/burden/
    """
    from hvantk.enrichex.pipeline import BurdenConfig, BurdenPipeline

    # Parse gene set collections (name:path pairs)
    collections = {}
    for pair in gene_sets_pairs:
        if ":" not in pair:
            click.echo(
                f"Error: --gene-sets must be 'name:path' format, got '{pair}'",
                err=True,
            )
            raise SystemExit(1)
        name, path = pair.split(":", 1)
        collections[name.strip()] = path.strip()

    # Parse variant classes
    vc_dict = {}
    if variant_classes:
        from hvantk.enrichex.burden import VariantFilter, build_variant_classes_from_presets

        class_names = [c.strip() for c in variant_classes.split(",")]
        base_filter = VariantFilter(
            max_af=max_af,
            min_score=None,
            consequences=None,
            pass_only=pass_only,
            min_gq=0,
            min_dp=0,
            af_field=af_field,
            score_field=score_field,
            consequence_field=consequence_field,
        )
        vc_dict = build_variant_classes_from_presets(class_names, base_filter=base_filter)

    covar_list = [c.strip() for c in covariates.split(",")] if covariates else []

    config = BurdenConfig(
        cohort_mt_path=cohort_mt,
        phenotype_ht_path=phenotypes or "",
        phenotype_field=phenotype_field,
        phenotype_type=phenotype_type,
        covariate_fields=covar_list,
        sample_id_field=sample_id_field,
        gene_set_collections=collections,
        variant_classes=vc_dict,
        gene_field=gene_field,
        genotype_aggregation=genotype_aggregation,
        normalize_by_length=normalize_by_length,
        gene_lengths_path=gene_lengths,
        min_carriers=min_carriers,
        correction_method=correction,
        alpha=alpha,
        competitive=competitive,
        n_permutations=n_permutations,
        permutation_seed=permutation_seed,
        output_dir=output_dir,
        generate_report=generate_report,
    )

    pipeline = BurdenPipeline(config)

    if dry_run:
        pipeline.show_plan()
        return

    # Initialize Hail before running
    from hvantk.core.hail_context import init_hail

    init_hail()

    combined_df = pipeline.run()

    click.echo(f"\nTotal results: {len(combined_df)} rows")
    click.echo(f"Output directory: {output_dir}")
