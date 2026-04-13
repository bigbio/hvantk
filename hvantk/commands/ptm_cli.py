"""
PTM CLI Commands - Post-Translational Modification Analysis

This module provides CLI commands for the PTM variant classification pipeline:
- build: Download PTM data, map coordinates, build Hail Table
- annotate: Annotate a variant table with PTM site information
- landscape: PTM-variant overlap analysis (Q1)
- export-strata: Export PTM/non-PTM variant strata for predictor evaluation (Q2)
- population: Population-level PTM-variant burden (Q3)
- report: Generate summary report
"""

import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


@click.group(
    name="ptm",
    help="Post-translational modification variant classification commands.",
    context_settings=CONTEXT_SETTINGS,
)
@click.pass_context
def ptm_group(ctx):
    """PTM command group for variant classification analysis.

    \b
    Workflow:
      1. hvantk ptm build          — Download PTM data + map coordinates + build Hail Table
      2. hvantk ptm annotate       — Annotate variants with PTM site information
      3. hvantk ptm landscape      — PTM-variant overlap analysis (Q1)
      4. hvantk ptm export-strata  — Export PTM/non-PTM variant lists for PSROC (Q2)
      5. hvantk ptm population     — Population-level PTM-variant burden (Q3)
      6. hvantk ptm constraint     — Stratified AF depletion at PTM codons (by tissue/cell-type)
      7. hvantk ptm report         — Generate summary report

    \b
    For predictor evaluation (Q2), compose with PSROC:
      hvantk ptm export-strata --annotated-ht clinvar_ptm.ht -o strata/
      hvantk psroc --variants strata/ptm_variants.txt --clinvar-ht clinvar.ht ...
      hvantk psroc --variants strata/non_ptm_variants.txt --clinvar-ht clinvar.ht ...

    \b
    Examples:
      hvantk ptm build --output-dir data/ptm/ --output-ht data/ptm/ptm_sites.ht
      hvantk ptm annotate --variants-ht clinvar.ht --ptm-ht ptm_sites.ht -o annotated.ht
    """
    ctx.ensure_object(dict)


@ptm_group.command("build")
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(),
    required=True,
    help="Directory for intermediate files (GTF, UniProt TSV, mapped TSV)",
)
@click.option(
    "--output-ht",
    type=str,
    required=True,
    help="Path to write the output PTM sites Hail Table (.ht)",
)
@click.option(
    "--gtf-path",
    type=click.Path(exists=True),
    default=None,
    help="Path to pre-downloaded Ensembl GTF (skips download if provided)",
)
@click.option(
    "--ptm-tsv",
    type=click.Path(exists=True),
    default=None,
    help="Path to pre-downloaded UniProt PTM TSV (skips API query if provided)",
)
@click.option(
    "--peptideatlas-tsv",
    type=click.Path(exists=True),
    default=None,
    help="Path to PeptideAtlas phospho intermediate TSV (from 'hvantk download peptideatlas-phospho')",
)
@click.option(
    "--cptac-tsv",
    type=click.Path(exists=True),
    default=None,
    help="Path to CPTAC phospho intermediate TSV (from 'hvantk download cptac-phospho')",
)
@click.option(
    "--flanking-codons",
    type=int,
    default=5,
    show_default=True,
    help="Number of flanking codons for proximal window",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing outputs",
)
@click.pass_context
def ptm_build(
    ctx, output_dir, output_ht, gtf_path, ptm_tsv, peptideatlas_tsv, cptac_tsv, flanking_codons, overwrite
):
    """Download PTM data, map coordinates to genome, and build a Hail Table.

    \b
    This is the main entry point for the PTM pipeline (Phases 1-2). It:
      1. Downloads the Ensembl GTF (if not provided via --gtf-path)
      2. Downloads UniProt PTM data (if not provided via --ptm-tsv)
      3. Parses the GTF and maps PTM sites to genomic coordinates
      4. Builds a Hail Table keyed by locus

    \b
    Examples:
      hvantk ptm build --output-dir data/ptm/ --output-ht data/ptm/ptm_sites.ht
      hvantk ptm build --gtf-path data/ref/Homo_sapiens.GRCh38.113.gtf.gz \\
                       --ptm-tsv data/ptm/uniprot-ptm-human.tsv \\
                       --output-ht data/ptm/ptm_sites.ht --output-dir data/ptm/
    """
    try:
        from hvantk.ptm.pipeline import PTMBuildConfig, ptm_build_pipeline

        config = PTMBuildConfig(
            output_dir=output_dir,
            output_ht=output_ht,
            gtf_path=gtf_path,
            ptm_tsv=ptm_tsv,
            peptideatlas_tsv=peptideatlas_tsv,
            cptac_tsv=cptac_tsv,
            flanking_codons=flanking_codons,
            overwrite=overwrite,
        )

        errors = config.validate()
        if errors:
            click.echo("Configuration validation failed:", err=True)
            for error in errors:
                click.echo(f"  - {error}", err=True)
            ctx.exit(1)

        result = ptm_build_pipeline(config)

        click.echo(
            f"Mapping complete: {result.n_mapped}/{result.n_total} mapped "
            f"({100 * result.n_mapped / max(result.n_total, 1):.1f}%), "
            f"{result.n_failed} failed"
        )
        click.echo(f"Resolution: {result.resolution_counts}")
        click.echo(f"Mapped TSV: {result.mapped_tsv_path}")
        click.echo(f"Hail Table: {result.output_ht}")

    except Exception as e:
        logger.exception(f"PTM build failed: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)


@ptm_group.command("annotate")
@click.option(
    "--variants-ht", type=str, required=True, help="Path to variant Hail Table"
)
@click.option("--ptm-ht", type=str, required=True, help="Path to PTM sites Hail Table")
@click.option(
    "-o", "--output-ht", type=str, required=True, help="Path to write annotated Table"
)
@click.option("--flanking-codons", type=int, default=5, show_default=True)
@click.option("--overwrite", is_flag=True)
@click.pass_context
def ptm_annotate(ctx, variants_ht, ptm_ht, output_ht, flanking_codons, overwrite):
    """Annotate a variant table with PTM site information.

    \b
    Cross-references variants with PTM sites using position-based joins.
    For each variant, determines whether it falls directly at a PTM codon
    or within the proximal flanking window.

    \b
    Output fields added to the variant table:
      is_ptm_site      — variant at a PTM-modified codon
      is_ptm_proximal  — variant within flanking window (not at codon)
      ptm_types        — set of PTM categories (e.g., phosphorylation)
      ptm_distance     — approximate distance in residues to nearest PTM site

    \b
    Examples:
      hvantk ptm annotate --variants-ht clinvar.ht --ptm-ht ptm_sites.ht -o annotated.ht
      hvantk ptm annotate --variants-ht gnomad.ht --ptm-ht ptm_sites.ht -o gnomad_ptm.ht --flanking-codons 7
    """
    try:
        import hail as hl
        from hvantk.ptm.annotate import annotate_variants_with_ptm

        variants = hl.read_table(variants_ht)
        ptm = hl.read_table(ptm_ht)

        result = annotate_variants_with_ptm(
            variants, ptm, flanking_codons=flanking_codons
        )
        result = result.checkpoint(output_ht, overwrite=overwrite)

        import hail as hl

        counts = result.aggregate(
            hl.struct(
                n_total=hl.agg.count(),
                n_ptm_site=hl.agg.filter(result.is_ptm_site, hl.agg.count()),
                n_proximal=hl.agg.filter(result.is_ptm_proximal, hl.agg.count()),
            )
        )

        click.echo(f"Annotated {counts.n_total:,} variants:")
        click.echo(f"  PTM site:  {counts.n_ptm_site:,}")
        click.echo(f"  Proximal:  {counts.n_proximal:,}")
        click.echo(f"Output: {output_ht}")

    except Exception as e:
        logger.exception(f"PTM annotation failed: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)


@ptm_group.command("landscape")
@click.option(
    "--clinvar-ht", type=str, required=True, help="Path to ClinVar Hail Table"
)
@click.option("--ptm-ht", type=str, required=True, help="Path to PTM sites Hail Table")
@click.option(
    "-o", "--output", type=click.Path(), required=True, help="Output directory"
)
@click.option("--flanking-codons", type=int, default=5, show_default=True)
@click.option("--save-plots", is_flag=True, help="Save plots alongside JSON output")
@click.pass_context
def ptm_landscape_cmd(ctx, clinvar_ht, ptm_ht, output, flanking_codons, save_plots):
    """PTM-variant overlap and enrichment analysis (Q1).

    \b
    Cross-references ClinVar P/LP variants with PTM sites, computes overlap
    counts per PTM category, enrichment (Fisher's exact test), and distance
    distribution.

    \b
    Output:
      landscape_summary.json — counts, enrichment, per-category overlaps
      *.png (with --save-plots) — landscape summary, category, distance plots

    \b
    Examples:
      hvantk ptm landscape --clinvar-ht clinvar.ht --ptm-ht ptm_sites.ht -o results/landscape/
      hvantk ptm landscape --clinvar-ht clinvar.ht --ptm-ht ptm_sites.ht -o results/ --save-plots
    """
    try:
        import os
        import hail as hl
        from hvantk.ptm.analysis import ptm_landscape

        clinvar = hl.read_table(clinvar_ht)
        ptm = hl.read_table(ptm_ht)

        result = ptm_landscape(clinvar, ptm, output, flanking_codons=flanking_codons)
        click.echo(result.summary())

        if save_plots:
            import matplotlib

            matplotlib.use("Agg")
            from hvantk.ptm.plot import (
                plot_landscape_summary,
                plot_overlap_by_category,
                plot_distance_distribution,
            )

            plot_landscape_summary(
                result, os.path.join(output, "landscape_summary.png")
            )
            plot_overlap_by_category(
                result, os.path.join(output, "overlap_by_category.png")
            )
            plot_distance_distribution(
                result, os.path.join(output, "distance_distribution.png")
            )
            click.echo(f"Plots saved to {output}")

    except Exception as e:
        logger.exception(f"PTM landscape failed: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)


@ptm_group.command("export-strata")
@click.option(
    "--annotated-ht",
    type=str,
    required=True,
    help="Path to PTM-annotated variant Hail Table",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    required=True,
    help="Output directory for variant lists",
)
@click.pass_context
def ptm_export_strata(ctx, annotated_ht, output):
    """Export PTM-stratified variant lists for downstream analysis (Q2).

    \b
    Takes a PTM-annotated variant table (output of 'hvantk ptm annotate') and
    exports two variant lists in chr:pos:ref:alt format:
      ptm_variants.txt     — variants at or proximal to PTM sites
      non_ptm_variants.txt — variants not near PTM sites

    \b
    These lists can be fed to hvantk psroc for stratified predictor evaluation:
      hvantk psroc --variants ptm_variants.txt --clinvar-ht clinvar.ht ...

    \b
    Examples:
      hvantk ptm annotate --variants-ht clinvar.ht --ptm-ht ptm_sites.ht -o clinvar_ptm.ht
      hvantk ptm export-strata --annotated-ht clinvar_ptm.ht -o strata/
    """
    try:
        import hail as hl
        from hvantk.ptm.analysis import export_ptm_strata

        ht = hl.read_table(annotated_ht)
        paths = export_ptm_strata(ht, output)

        for name, path in paths.items():
            click.echo(f"  {name}: {path}")

    except Exception as e:
        logger.exception(f"PTM export-strata failed: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)


@ptm_group.command("population")
@click.option("--gnomad-ht", type=str, required=True, help="Path to gnomAD Hail Table")
@click.option("--ptm-ht", type=str, required=True, help="Path to PTM sites Hail Table")
@click.option(
    "--ccr-ht", type=str, default=None, help="Path to CCR Hail Table (optional)"
)
@click.option(
    "--af-field",
    type=str,
    default="AF",
    show_default=True,
    help="AF field name in gnomAD table",
)
@click.option(
    "-o", "--output", type=click.Path(), required=True, help="Output directory"
)
@click.option("--flanking-codons", type=int, default=5, show_default=True)
@click.option(
    "--save-plots", is_flag=True, help="Save AF comparison plot alongside JSON output"
)
@click.pass_context
def ptm_population_cmd(
    ctx, gnomad_ht, ptm_ht, ccr_ht, af_field, output, flanking_codons, save_plots
):
    """Population-level PTM-variant allele frequency analysis (Q3).

    \b
    Compares allele frequency distributions at PTM sites vs non-PTM coding
    positions in gnomAD. Identifies PTM sites under purifying selection
    (zero or near-zero AF). Optionally cross-references with CCR scores.

    \b
    Output:
      population_summary.json — AF statistics, counts, CCR comparison
      *.png (with --save-plots) — AF comparison plot

    \b
    Examples:
      hvantk ptm population --gnomad-ht gnomad.ht --ptm-ht ptm_sites.ht -o results/population/
      hvantk ptm population --gnomad-ht gnomad.ht --ptm-ht ptm_sites.ht --ccr-ht ccr.ht -o results/population/
    """
    try:
        import os
        import hail as hl
        from hvantk.ptm.analysis import ptm_population

        gnomad = hl.read_table(gnomad_ht)
        ptm = hl.read_table(ptm_ht)
        ccr = hl.read_table(ccr_ht) if ccr_ht else None

        result = ptm_population(
            gnomad,
            ptm,
            output,
            ccr_ht=ccr,
            af_field=af_field,
            flanking_codons=flanking_codons,
        )
        click.echo(result.summary())

        if save_plots:
            import matplotlib

            matplotlib.use("Agg")
            from hvantk.ptm.plot import plot_population_af

            plot_population_af(result, os.path.join(output, "population_af.png"))
            click.echo(f"Plot saved to {output}")

    except Exception as e:
        logger.exception(f"PTM population analysis failed: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)


@ptm_group.command("report")
@click.option(
    "-o", "--output", type=click.Path(), required=True, help="Output HTML report path"
)
@click.option(
    "--landscape-json",
    type=click.Path(exists=True),
    default=None,
    help="Path to landscape_summary.json (from 'hvantk ptm landscape')",
)
@click.option(
    "--population-json",
    type=click.Path(exists=True),
    default=None,
    help="Path to population_summary.json (from 'hvantk ptm population')",
)
@click.option("--title", type=str, default="PTM-Variant Analysis Report")
@click.option("--description", type=str, default=None)
@click.pass_context
def ptm_report(ctx, output, landscape_json, population_json, title, description):
    """Generate PTM analysis summary HTML report.

    \b
    Combines landscape (Q1) and/or population (Q3) results into a single
    HTML report with embedded plots and summary tables.

    \b
    Examples:
      hvantk ptm report -o report.html --landscape-json results/landscape/landscape_summary.json
      hvantk ptm report -o report.html \\
          --landscape-json results/landscape/landscape_summary.json \\
          --population-json results/population/population_summary.json
    """
    if not landscape_json and not population_json:
        click.echo(
            "Error: at least one of --landscape-json or --population-json is required.",
            err=True,
        )
        ctx.exit(1)

    try:
        import json
        from hvantk.ptm.analysis import PTMLandscapeResult, PTMPopulationResult
        from hvantk.ptm.report import generate_report

        landscape_result = None
        if landscape_json:
            with open(landscape_json) as f:
                data = json.load(f)
            landscape_result = PTMLandscapeResult(
                n_variants=data.get("n_variants", 0),
                n_pathogenic=data.get("n_pathogenic", 0),
                n_benign=data.get("n_benign", 0),
                n_ptm_site_pathogenic=data.get("ptm_site", {}).get("pathogenic", 0),
                n_ptm_proximal_pathogenic=data.get("ptm_proximal", {}).get(
                    "pathogenic", 0
                ),
                n_ptm_site_benign=data.get("ptm_site", {}).get("benign", 0),
                n_ptm_proximal_benign=data.get("ptm_proximal", {}).get("benign", 0),
                enrichment_odds_ratio=data.get("enrichment", {}).get("odds_ratio", 0.0),
                enrichment_ci_low=data.get("enrichment", {}).get("ci_low", 0.0),
                enrichment_ci_high=data.get("enrichment", {}).get(
                    "ci_high", float("inf")
                ),
                enrichment_p_value=data.get("enrichment", {}).get("p_value", 1.0),
                overlap_by_category=data.get("overlap_by_category", {}),
                category_enrichment=data.get("category_enrichment", {}),
                distance_distribution={
                    int(k): v for k, v in data.get("distance_distribution", {}).items()
                },
            )

        population_result = None
        if population_json:
            with open(population_json) as f:
                data = json.load(f)
            population_result = PTMPopulationResult(
                n_variants=data.get("n_variants", 0),
                n_ptm_site=data.get("n_ptm_site", 0),
                n_ptm_proximal=data.get("n_ptm_proximal", 0),
                n_non_ptm=data.get("n_non_ptm", 0),
                mean_af_ptm_site=data.get("mean_af", {}).get("ptm_site", 0.0),
                mean_af_ptm_proximal=data.get("mean_af", {}).get("ptm_proximal", 0.0),
                mean_af_non_ptm=data.get("mean_af", {}).get("non_ptm", 0.0),
                n_zero_af_ptm=data.get("n_zero_af_ptm", 0),
                ptm_site_afs=data.get("ptm_site_afs", []),
                proximal_afs=data.get("proximal_afs", []),
                non_ptm_afs=data.get("non_ptm_afs", []),
                ccr_mean_ptm=data.get("ccr_mean_ptm"),
                ccr_mean_non_ptm=data.get("ccr_mean_non_ptm"),
            )

        generate_report(
            output_path=output,
            landscape_result=landscape_result,
            population_result=population_result,
            title=title,
            description=description,
        )
        click.echo(f"Report: {output}")

    except Exception as e:
        logger.exception(f"PTM report failed: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)


@ptm_group.command("constraint")
@click.option(
    "--variants-ht",
    type=str,
    required=True,
    help="Path to PTM-annotated variant Hail Table (from `hvantk ptm annotate`).",
)
@click.option(
    "--expression-source",
    type=click.Choice(["hail-mt", "anndata", "tabular"], case_sensitive=False),
    required=True,
    help="Expression data backend.",
)
@click.option(
    "--expression-path",
    type=click.Path(exists=True),
    required=True,
    help="Path to the expression data (.mt / .h5ad / .parquet|.pkl|.tsv).",
)
@click.option(
    "--grouping",
    type=str,
    required=True,
    help="Metadata field to stratify by (e.g. SMTSD, major_cell_class).",
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    required=True,
    help="Directory for TSVs, plots, and HTML report.",
)
@click.option(
    "--label-filter",
    type=click.Choice(["TN", "TP", "all"], case_sensitive=False),
    default="TN",
    show_default=True,
    help="Label subset to include (default: TN = benign/likely-benign).",
)
@click.option("--label-field", type=str, default="rf_label", show_default=True)
@click.option("--gene-field", type=str, default="gene_symbol", show_default=True)
@click.option("--af-field", type=str, default="gnomad_af_genomes", show_default=True)
@click.option("--loeuf-field", type=str, default="loeuf", show_default=True)
@click.option("--ptm-category-field", type=str, default="ptm_types", show_default=True)
@click.option(
    "--gene-id-mapping",
    type=click.Path(exists=True),
    default=None,
    help="Optional two-column TSV mapping expression gene IDs to symbols.",
)
@click.option(
    "--expression-metric",
    type=click.Choice(["median", "mean", "median_nonzero"], case_sensitive=False),
    default="median",
    show_default=True,
)
@click.option(
    "--min-cells-per-group",
    type=int,
    default=50,
    show_default=True,
    help="Drop groups with fewer cells/samples (AnnData/Hail MT only).",
)
@click.option(
    "--min-variants-per-group",
    type=int,
    default=20,
    show_default=True,
    help="Skip groups with fewer than N variants across PTM + non-PTM.",
)
@click.option("--flanking-codons", type=int, default=7, show_default=True)
@click.option(
    "--expressed-threshold",
    type=float,
    default=1.0,
    show_default=True,
    help="Keep only genes whose max group expression meets this floor.",
)
@click.option("--overwrite", is_flag=True)
@click.pass_context
def ptm_constraint(
    ctx,
    variants_ht,
    expression_source,
    expression_path,
    grouping,
    output_dir,
    label_filter,
    label_field,
    gene_field,
    af_field,
    loeuf_field,
    ptm_category_field,
    gene_id_mapping,
    expression_metric,
    min_cells_per_group,
    min_variants_per_group,
    flanking_codons,
    expressed_threshold,
    overwrite,
):
    """Stratified PTM constraint analysis across groups (tissue / cell type).

    \b
    Compares gnomAD allele-frequency distributions between PTM-proximal and
    non-PTM variants, stratified by a metadata field from the expression
    dataset. Runs five tests (per-group ranking, tau quartile, LOEUF x group
    factorial, PTM category x group heatmap, within-gene Wilcoxon) and emits
    TSVs, four PNG panels, and one self-contained HTML report.

    \b
    Note: this is a stratified DEPLETION analysis — not a per-variant scorer.
    For per-site PTM flags use `hvantk ptm annotate`.

    \b
    Examples:
      hvantk ptm constraint \\
          --variants-ht clinvar_ptm.ht \\
          --expression-source hail-mt \\
          --expression-path /data/GTEX_v7.mt \\
          --grouping SMTSD \\
          --output-dir results/ptm-gtex/

      hvantk ptm constraint \\
          --variants-ht clinvar_ptm.ht \\
          --expression-source anndata \\
          --expression-path /data/farah_2024.h5ad \\
          --grouping major_cell_class \\
          --output-dir results/ptm-farah/
    """
    try:
        from hvantk.ptm.constraint import PTMConstraintConfig, run_ptm_constraint

        config = PTMConstraintConfig(
            variants_ht_path=variants_ht,
            expression_source=expression_source.lower(),
            expression_path=expression_path,
            grouping=grouping,
            output_dir=output_dir,
            label_filter=label_filter.upper() if label_filter != "all" else "all",
            label_field=label_field,
            gene_field=gene_field,
            af_field=af_field,
            loeuf_field=loeuf_field,
            ptm_category_field=ptm_category_field,
            gene_id_mapping=gene_id_mapping,
            expression_metric=expression_metric.lower(),
            min_cells_per_group=min_cells_per_group,
            min_variants_per_group=min_variants_per_group,
            flanking_codons=flanking_codons,
            expressed_threshold=expressed_threshold,
            overwrite=overwrite,
        )

        errors = config.validate()
        if errors:
            click.echo("Configuration validation failed:", err=True)
            for err in errors:
                click.echo(f"  - {err}", err=True)
            ctx.exit(1)

        result = run_ptm_constraint(config)
        click.echo(result.summary())
        click.echo(f"Output: {config.output_dir}")

    except Exception as e:
        logger.exception(f"PTM constraint failed: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
