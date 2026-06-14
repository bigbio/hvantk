"""
CLI commands for QTL cascade analysis.

Commands::

    hvantk qtlcascade cascade  — Build the eQTL ⊕ pQTL cascade join
    hvantk qtlcascade coloc    — Run colocalization on cascade genes
    hvantk qtlcascade run      — Full pipeline (cascade + coloc + gene summary + report)
    hvantk qtlcascade report   — Generate HTML report from existing results
"""

import logging
import math

import click

from hvantk.core.config import CONTEXT_SETTINGS
from hvantk.algorithms.qtlcascade.constants import (
    DEFAULT_COLOC_H4_THRESHOLD,
    DEFAULT_COLOC_WINDOW_KB,
    EQTL_CATALOGUE_DEFAULT_STUDY,
    DEFAULT_FINEMAP_SUPERPOP,
)

logger = logging.getLogger(__name__)


@click.group(
    name="qtlcascade",
    help="Molecular QTL cascade analysis (eQTL → pQTL → disease).",
    context_settings=CONTEXT_SETTINGS,
)
@click.pass_context
def qtlcascade_group(ctx):
    """QTL cascade command group."""
    ctx.ensure_object(dict)


# ---------------------------------------------------------------------------
# cascade — build the join
# ---------------------------------------------------------------------------


@qtlcascade_group.command("cascade")
@click.option("--eqtl-ht", required=True, type=str, help="Path to eQTL Hail Table.")
@click.option("--pqtl-ht", required=True, type=str, help="Path to pQTL Hail Table.")
@click.option(
    "-o", "--output", required=True, type=str, help="Output cascade Hail Table path."
)
@click.option("--tissue", type=str, default=None, help="Filter to this tissue.")
@click.option(
    "--eqtl-p",
    type=float,
    default=5e-8,
    show_default=True,
    help="eQTL p-value threshold.",
)
@click.option(
    "--pqtl-p",
    type=float,
    default=5e-8,
    show_default=True,
    help="pQTL p-value threshold.",
)
@click.option("--overwrite", is_flag=True, help="Overwrite existing output.")
@click.pass_context
def cascade_cmd(ctx, eqtl_ht, pqtl_ht, output, tissue, eqtl_p, pqtl_p, overwrite):
    """Build the QTL cascade (eQTL ⊕ pQTL outer join + classification)."""
    from hvantk.core.utils.hail_context import init_hail

    init_hail()

    from hvantk.algorithms.qtlcascade.cascade import build_cascade

    ht = build_cascade(
        eqtl_ht_path=eqtl_ht,
        pqtl_ht_path=pqtl_ht,
        output_path=output,
        eqtl_p_threshold=eqtl_p,
        pqtl_p_threshold=pqtl_p,
        tissue=tissue,
        overwrite=overwrite,
    )
    click.echo(f"Cascade table written to {output}")
    ht.describe()


# ---------------------------------------------------------------------------
# coloc — colocalization
# ---------------------------------------------------------------------------


@qtlcascade_group.command("coloc")
@click.option(
    "--eqtl-allpairs", required=True, type=str, help="Path to allpairs eQTL Hail Table."
)
@click.option(
    "--pqtl-allpairs", required=True, type=str, help="Path to allpairs pQTL Hail Table."
)
@click.option(
    "--cascade-genes",
    required=True,
    type=str,
    help="File with one gene_id per line (cascade genes).",
)
@click.option(
    "--tissue", type=str, default=None, help="Filter allpairs to this tissue."
)
@click.option(
    "--window-kb",
    type=int,
    default=500,
    show_default=True,
    help="Regional window (±kb) for coloc.",
)
@click.option(
    "-o", "--output", required=True, type=str, help="Output TSV path for coloc results."
)
@click.pass_context
def coloc_cmd(
    ctx, eqtl_allpairs, pqtl_allpairs, cascade_genes, tissue, window_kb, output
):
    """Run colocalization ABF on cascade genes."""
    from pathlib import Path

    from hvantk.core.utils.hail_context import init_hail

    init_hail()

    from hvantk.algorithms.qtlcascade.coloc import run_coloc_per_gene

    genes = Path(cascade_genes).read_text().strip().splitlines()
    genes = [g.strip() for g in genes if g.strip()]
    click.echo(f"Running coloc for {len(genes)} genes")

    df = run_coloc_per_gene(
        eqtl_allpairs_ht_path=eqtl_allpairs,
        pqtl_allpairs_ht_path=pqtl_allpairs,
        cascade_genes=genes,
        tissue=tissue,
        window_kb=window_kb,
    )

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, sep="\t", index=False)
    click.echo(f"Coloc results written to {output} ({len(df)} genes)")


# ---------------------------------------------------------------------------
# gwas-coloc — GWAS → effector coloc (FinnGen × eQTL Catalogue) + SuSiE confirm
# ---------------------------------------------------------------------------


@qtlcascade_group.command("gwas-coloc")
@click.option("--endpoint", required=True, type=str,
              help="FinnGen R10 endpoint code (e.g. I9_AF).")
@click.option("--chrom", required=True, type=str, help="Chromosome (GRCh38, no 'chr').")
@click.option("--lead", required=True, type=int, help="Lead variant position (GRCh38).")
@click.option("--eqtl", "eqtl_dataset", required=True, type=str,
              help="eQTL Catalogue dataset id / URL / local tabix (e.g. QTD000251).")
@click.option("--eqtl-study", type=str, default=EQTL_CATALOGUE_DEFAULT_STUDY,
              show_default=True, help="eQTL Catalogue study id.")
@click.option("--window-kb", type=int, default=DEFAULT_COLOC_WINDOW_KB,
              show_default=True, help="Regional window (±kb) around the lead.")
@click.option("--gene", "gene_of_interest", type=str, default=None,
              help="ENSG of interest to confirm (default: the ABF-top gene).")
@click.option("--fine-map/--no-fine-map", default=True, show_default=True,
              help="Run SuSiE/coloc.susie confirmation (needs R+susieR+coloc, bcftools, curl).")
@click.option("--gwas-n", "gwas_n", type=int, default=None,
              help="GWAS sample size (required for fine-mapping).")
@click.option("--eqtl-n", "eqtl_n", type=int, default=None,
              help="eQTL sample size (required for fine-mapping).")
@click.option("--superpop", type=str, default=DEFAULT_FINEMAP_SUPERPOP,
              show_default=True, help="1000G super-population for the LD reference.")
@click.option("--ld-cache-dir", type=str, default=None,
              help="Cache directory for the 1000G LD reference.")
@click.option("-o", "--output-dir", required=True, type=str, help="Output directory.")
@click.pass_context
def gwas_coloc_cmd(ctx, endpoint, chrom, lead, eqtl_dataset, eqtl_study, window_kb,
                   gene_of_interest, fine_map, gwas_n, eqtl_n, superpop,
                   ld_cache_dir, output_dir):
    """GWAS → effector colocalization with optional SuSiE fine-map confirmation.

    Ranks cis effectors at a GWAS locus by ABF H4, then (default) confirms the
    lead effector with SuSiE/coloc.susie to separate genuine colocalization
    (CONFIRMED) from single-variant-ABF artifacts (REFUTED).
    """
    from hvantk.algorithms.qtlcascade.gwas_pipeline import (
        GwasColocConfig,
        run_gwas_coloc_pipeline,
    )

    config = GwasColocConfig(
        endpoint=endpoint, chrom=chrom, lead=lead, eqtl_dataset=eqtl_dataset,
        eqtl_study=eqtl_study, window_kb=window_kb, gene_of_interest=gene_of_interest,
        fine_map=fine_map, gwas_N=gwas_n, eqtl_N=eqtl_n, superpop=superpop,
        ld_cache_dir=ld_cache_dir, output_dir=output_dir,
    )
    errors = config.validate()
    if errors:
        for e in errors:
            click.echo(f"  ERROR: {e}", err=True)
        ctx.exit(1)

    report = run_gwas_coloc_pipeline(config)
    res = report["results"]
    gmp = report["gwas"]["min_p_in_region"]
    gmp_str = f"{gmp:.1e}" if isinstance(gmp, float) and math.isfinite(gmp) else "n/a"
    click.echo(f"\nRegion {report['region']}  |  GWAS min-p {gmp_str}  |  "
               f"{res['n_genes_tested']} genes tested")
    click.echo(f"Top effector: {res['top_effector']} (PP4={res['top_PP4']})")
    if report["fine_map"]:
        fmr = report["fine_map"]
        if fmr["available"]:
            click.echo(f"Fine-map: credible sets GWAS={fmr['credible_sets_gwas']} "
                       f"eQTL={fmr['credible_sets_eqtl']}; "
                       f"coloc.susie PP4={fmr['coloc_susie_PP4']}")
        else:
            click.echo(f"Fine-map: {fmr['note']}", err=True)
    click.echo(f"VERDICT: {report['verdict']}")
    click.echo(f"Report: {res['report_json']}")


# ---------------------------------------------------------------------------
# run — full pipeline
# ---------------------------------------------------------------------------


@qtlcascade_group.command("run")
@click.option("--eqtl-ht", required=True, type=str, help="Path to eQTL Hail Table.")
@click.option("--pqtl-ht", required=True, type=str, help="Path to pQTL Hail Table.")
@click.option(
    "--eqtl-allpairs",
    type=str,
    default="",
    help="Path to allpairs eQTL HT (for coloc).",
)
@click.option(
    "--pqtl-allpairs",
    type=str,
    default="",
    help="Path to allpairs pQTL HT (for coloc).",
)
@click.option(
    "--constraint-ht",
    type=str,
    default="",
    help="Path to gnomAD constraint HT (LOEUF overlay).",
)
@click.option(
    "--disease-genes-ht", type=str, default="", help="Path to disease-gene HT."
)
@click.option(
    "--tissues",
    type=str,
    default="",
    help="Comma-separated tissue list for multi-tissue run.",
)
@click.option("-o", "--output-dir", required=True, type=str, help="Output directory.")
@click.option(
    "--eqtl-p",
    type=float,
    default=5e-8,
    show_default=True,
    help="eQTL p-value threshold.",
)
@click.option(
    "--pqtl-p",
    type=float,
    default=5e-8,
    show_default=True,
    help="pQTL p-value threshold.",
)
@click.option(
    "--window-kb", type=int, default=500, show_default=True, help="Coloc window (±kb)."
)
@click.option("--no-plots", is_flag=True, help="Skip plot generation.")
@click.option("--no-report", is_flag=True, help="Skip HTML report.")
@click.option("--overwrite", is_flag=True, help="Overwrite existing outputs.")
@click.option("--dry-run", is_flag=True, help="Show plan without executing.")
@click.pass_context
def run_cmd(
    ctx,
    eqtl_ht,
    pqtl_ht,
    eqtl_allpairs,
    pqtl_allpairs,
    constraint_ht,
    disease_genes_ht,
    tissues,
    output_dir,
    eqtl_p,
    pqtl_p,
    window_kb,
    no_plots,
    no_report,
    overwrite,
    dry_run,
):
    """Run the full QTL cascade pipeline."""
    from hvantk.algorithms.qtlcascade.pipeline import CascadeConfig, CascadePipeline

    tissue_list = [t.strip() for t in tissues.split(",") if t.strip()]

    config = CascadeConfig(
        eqtl_ht=eqtl_ht,
        pqtl_ht=pqtl_ht,
        eqtl_allpairs_ht=eqtl_allpairs,
        pqtl_allpairs_ht=pqtl_allpairs,
        constraint_ht=constraint_ht,
        disease_genes_ht=disease_genes_ht,
        tissues=tissue_list,
        output_dir=output_dir,
        eqtl_p_threshold=eqtl_p,
        pqtl_p_threshold=pqtl_p,
        coloc_window_kb=window_kb,
        generate_plots=not no_plots,
        generate_report=not no_report,
        overwrite=overwrite,
    )

    errors = config.validate()
    if errors:
        for e in errors:
            click.echo(f"  ERROR: {e}", err=True)
        ctx.exit(1)

    pipeline = CascadePipeline(config)

    if dry_run:
        plan = pipeline.show_plan()
        click.echo(plan)
        return

    click.echo("\n" + "=" * 70)
    click.echo("Starting QTL Cascade Pipeline")
    click.echo("=" * 70 + "\n")

    if tissue_list:
        results = pipeline.run_collection()
        click.echo(f"\nCompleted: {len(results)} tissue(s)")
        for tissue, res in results.items():
            n_pass = 0
            if res.coloc_df is not None and not res.coloc_df.empty:
                n_pass = (res.coloc_df["H4"] > DEFAULT_COLOC_H4_THRESHOLD).sum()
            click.echo(
                f"  {tissue}: {res.n_cascade_genes} genes, " f"{n_pass} colocalised"
            )
    else:
        result = pipeline.run()
        click.echo(f"\nCascade genes: {result.n_cascade_genes}")
        click.echo(f"Class counts: {result.class_counts}")

    click.echo(f"\nOutputs: {output_dir}")


# ---------------------------------------------------------------------------
# report — generate from existing results
# ---------------------------------------------------------------------------


@qtlcascade_group.command("report")
@click.option("--gene-summary", type=str, default=None, help="Gene summary TSV.")
@click.option("--coloc-results", type=str, default=None, help="Coloc results TSV.")
@click.option(
    "--plots-dir", type=str, default=None, help="Directory containing plot PNGs."
)
@click.option("-o", "--output", required=True, type=str, help="Output HTML path.")
@click.option(
    "--title", type=str, default="QTL Cascade Analysis Report", help="Report title."
)
@click.pass_context
def report_cmd(ctx, gene_summary, coloc_results, plots_dir, output, title):
    """Generate an HTML report from existing cascade results."""
    import pandas as pd
    from pathlib import Path
    from hvantk.algorithms.qtlcascade.report import generate_report

    gene_df = pd.read_csv(gene_summary, sep="\t") if gene_summary else None
    coloc_df = pd.read_csv(coloc_results, sep="\t") if coloc_results else None

    plot_paths = {}
    if plots_dir:
        for p in Path(plots_dir).glob("*.png"):
            plot_paths[p.stem] = str(p)

    generate_report(
        output_path=output,
        gene_summary_df=gene_df,
        coloc_df=coloc_df,
        plot_paths=plot_paths,
        title=title,
    )
    click.echo(f"Report written to {output}")
