"""
PTM CLI Commands - Post-Translational Modification Analysis

This module provides CLI commands for the PTM variant classification pipeline:
- build: Download PTM data, map coordinates, build Hail Table
- annotate: Annotate a variant table with PTM site information
- landscape: PTM-variant overlap analysis (Q1)
- evaluate: Predictor performance at PTM sites (Q2)
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
      1. hvantk ptm build       — Download PTM data + map coordinates + build Hail Table
      2. hvantk ptm annotate    — Annotate variants with PTM site information
      3. hvantk ptm landscape   — PTM-variant overlap analysis (Q1)
      4. hvantk ptm evaluate    — Predictor performance at PTM sites (Q2)
      5. hvantk ptm population  — Population-level PTM-variant burden (Q3)
      6. hvantk ptm report      — Generate summary report

    \b
    Examples:
      hvantk ptm build --output-dir data/ptm/ --output-ht data/ptm/ptm_sites.ht
      hvantk ptm annotate --variants-ht clinvar.ht --ptm-ht ptm_sites.ht -o annotated.ht
    """
    ctx.ensure_object(dict)


@ptm_group.command("build")
@click.option(
    "--output-dir", "-o",
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
    "--flanking-codons",
    type=int,
    default=5,
    show_default=True,
    help="Number of flanking codons for proximal window",
)
@click.option(
    "--overwrite", is_flag=True, help="Overwrite existing outputs",
)
@click.pass_context
def ptm_build(ctx, output_dir, output_ht, gtf_path, ptm_tsv, flanking_codons, overwrite):
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
            flanking_codons=flanking_codons,
            overwrite=overwrite,
        )

        errors = config.validate()
        if errors:
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
@click.option("--variants-ht", type=str, required=True, help="Path to variant Hail Table")
@click.option("--ptm-ht", type=str, required=True, help="Path to PTM sites Hail Table")
@click.option("-o", "--output-ht", type=str, required=True, help="Path to write annotated Table")
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

        n_total = result.count()
        n_ptm_site = result.filter(result.is_ptm_site).count()
        n_proximal = result.filter(result.is_ptm_proximal).count()

        click.echo(f"Annotated {n_total:,} variants:")
        click.echo(f"  PTM site:  {n_ptm_site:,}")
        click.echo(f"  Proximal:  {n_proximal:,}")
        click.echo(f"Output: {output_ht}")

    except Exception as e:
        logger.exception(f"PTM annotation failed: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)


@ptm_group.command("landscape")
@click.option("--clinvar-ht", type=str, required=True)
@click.option("--ptm-ht", type=str, required=True)
@click.option("-o", "--output", type=click.Path(), required=True)
@click.pass_context
def ptm_landscape(ctx, clinvar_ht, ptm_ht, output):
    """PTM-variant overlap analysis (Q1)."""
    click.echo("PTM landscape analysis not yet implemented (Phase 4a).")
    ctx.exit(0)


@ptm_group.command("evaluate")
@click.option("--clinvar-ht", type=str, required=True)
@click.option("--dbnsfp-ht", type=str, required=True)
@click.option("--ptm-ht", type=str, required=True)
@click.option("-o", "--output", type=click.Path(), required=True)
@click.pass_context
def ptm_evaluate(ctx, clinvar_ht, dbnsfp_ht, ptm_ht, output):
    """Predictor performance at PTM sites (Q2)."""
    click.echo("PTM predictor evaluation not yet implemented (Phase 4b).")
    ctx.exit(0)


@ptm_group.command("population")
@click.option("--gnomad-ht", type=str, required=True)
@click.option("--ptm-ht", type=str, required=True)
@click.option("--ccr-ht", type=str, default=None)
@click.option("-o", "--output", type=click.Path(), required=True)
@click.pass_context
def ptm_population(ctx, gnomad_ht, ptm_ht, ccr_ht, output):
    """Population-level PTM-variant burden (Q3)."""
    click.echo("PTM population analysis not yet implemented (Phase 4c).")
    ctx.exit(0)


@ptm_group.command("report")
@click.option("-o", "--output", type=click.Path(), required=True)
@click.pass_context
def ptm_report(ctx, output):
    """Generate PTM analysis summary report."""
    click.echo("PTM report generation not yet implemented (Phase 5).")
    ctx.exit(0)
