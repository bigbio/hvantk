"""
Flexible CLI to build individual Hail Tables/MatrixTables from raw inputs.

Examples:
  hvantk mktable clinvar --raw-input /path/to/clinvar.vcf.bgz --output-ht /out/clinvar.ht --ref-genome GRCh38
  hvantk mktable interactome --raw-input /path/to/interactome.bed.bgz --output-ht /out/interactome.ht
  hvantk mktable gevir --raw-input /path/to/gevir.tsv.bgz --output-ht /out/gevir.ht --fields oe_syn_upper,oe_mis_upper
  hvantk mktable gnomad-metrics --raw-input /path/to/gnomad.tsv.bgz --output-ht /out/gnomad.ht
  hvantk mktable ensembl-gene --raw-input /path/to/biomart.tsv.bgz --output-ht /out/ensembl.ht --no-canonical
  hvantk mktable dbnsfp --raw-input /path/to/dbNSFP4_variant.bgz --output-ht /out/dbnsfp.ht
"""

import logging
from typing import Optional, List

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


# Internal wrappers so tests can mock without importing hail-heavy modules at import time


def _create_clinvar_tb(*args, **kwargs):
    from hvantk.tables.table_builders import create_clinvar_tb

    return create_clinvar_tb(*args, **kwargs)


def _create_interactome_tb(*args, **kwargs):
    from hvantk.tables.table_builders import create_interactome_tb

    return create_interactome_tb(*args, **kwargs)


def _create_gevir_tb(*args, **kwargs):
    from hvantk.tables.table_builders import create_gevir_tb

    return create_gevir_tb(*args, **kwargs)


def _create_gnomad_constraint_gene_metrics_tb(*args, **kwargs):
    from hvantk.tables.table_builders import create_gnomad_constraint_gene_metrics_tb

    return create_gnomad_constraint_gene_metrics_tb(*args, **kwargs)


def _create_ensembl_gene_tb(*args, **kwargs):
    from hvantk.tables.table_builders import create_ensembl_gene_tb

    return create_ensembl_gene_tb(*args, **kwargs)


def _create_dbnsfp_tb(*args, **kwargs):
    from hvantk.tables.table_builders import create_dbnsfp_tb

    return create_dbnsfp_tb(*args, **kwargs)


def _create_clingen_gene_disease_tb(*args, **kwargs):
    from hvantk.tables.table_builders import create_clingen_gene_disease_tb

    return create_clingen_gene_disease_tb(*args, **kwargs)


def _create_hgnc_gene_tb(*args, **kwargs):
    from hvantk.tables.table_builders import create_hgnc_gene_tb

    return create_hgnc_gene_tb(*args, **kwargs)


@click.group("mktable", context_settings=CONTEXT_SETTINGS)
def mktable_group():
    """Create a single annotation Table/MatrixTable from a raw input file."""
    pass


# Shared options
_raw_input_opt = click.option(
    "--raw-input",
    required=True,
    type=str,
    help="Path to the raw input file (VCF/TSV/BED/etc.)",
)
_output_ht_opt = click.option(
    "--output-ht",
    required=True,
    type=str,
    help="Path to write the output Hail Table (.ht)",
)
_overwrite_opt = click.option(
    "--overwrite", is_flag=True, help="Overwrite existing outputs if present"
)
_export_tsv_opt = click.option(
    "--export-tsv",
    is_flag=True,
    help="Additionally export a flattened TSV (.tsv.bgz) next to the HT",
)
_ref_genome_opt = click.option(
    "--ref-genome",
    type=click.Choice(["GRCh38", "GRCh37"], case_sensitive=False),
    default="GRCh38",
    show_default=True,
    help="Reference genome to use for interval/locus parsing",
)


@mktable_group.command("clinvar")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
@_ref_genome_opt
def mktable_clinvar(
    raw_input: str, output_ht: str, overwrite: bool, export_tsv: bool, ref_genome: str
):
    """Build a ClinVar Hail Table from a VCF (keyed by locus, alleles)."""
    logger.info("Building ClinVar table")
    _create_clinvar_tb(
        input_path=raw_input,
        output_path=output_ht,
        overwrite=overwrite,
        export_tsv=export_tsv,
        reference_genome=ref_genome,
    )
    click.echo(f"ClinVar table created at {output_ht}")


@mktable_group.command("interactome")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
@_ref_genome_opt
def mktable_interactome(
    raw_input: str, output_ht: str, overwrite: bool, export_tsv: bool, ref_genome: str
):
    """Build an interactome Table from a BED (keyed by interval)."""
    logger.info("Building interactome table")
    _create_interactome_tb(
        input_path=raw_input,
        output_path=output_ht,
        overwrite=overwrite,
        export_tsv=export_tsv,
        reference_genome=ref_genome,
    )
    click.echo(f"Interactome table created at {output_ht}")


@mktable_group.command("gevir")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
@click.option(
    "--fields",
    type=str,
    default=None,
    help="Comma-separated list of fields to retain (optional)",
)
def mktable_gevir(
    raw_input: str,
    output_ht: str,
    overwrite: bool,
    export_tsv: bool,
    fields: Optional[str],
):
    """Build a GeVIR gene-level Table from a TSV (keyed by gene_id)."""
    selected: Optional[List[str]] = (
        [f.strip() for f in fields.split(",")] if fields else None
    )
    logger.info("Building GeVIR table")
    _create_gevir_tb(
        input_path=raw_input,
        output_path=output_ht,
        fields=selected,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"GEVIR table created at {output_ht}")


@mktable_group.command("gnomad-metrics")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
@click.option(
    "--fields",
    type=str,
    default=None,
    help="Comma-separated list of fields to retain (optional)",
)
def mktable_gnomad_metrics(
    raw_input: str,
    output_ht: str,
    overwrite: bool,
    export_tsv: bool,
    fields: Optional[str],
):
    """Build a gnomAD constraint metrics gene Table from a TSV (keyed by gene_id)."""
    selected: Optional[List[str]] = (
        [f.strip() for f in fields.split(",")] if fields else None
    )
    logger.info("Building gnomAD metrics table")
    _create_gnomad_constraint_gene_metrics_tb(
        input_path=raw_input,
        output_path=output_ht,
        fields=selected,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"gnomAD metrics table created at {output_ht}")


@mktable_group.command("ensembl-gene")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
@click.option(
    "--fields",
    type=str,
    default=None,
    help="Comma-separated list of fields to retain (optional)",
)
@click.option(
    "--canonical/--no-canonical",
    default=True,
    show_default=True,
    help="Include only canonical transcripts",
)
def mktable_ensembl_gene(
    raw_input: str,
    output_ht: str,
    overwrite: bool,
    export_tsv: bool,
    fields: Optional[str],
    canonical: bool,
):
    """Build an Ensembl gene annotation Table from a Biomart TSV (keyed by gene_id)."""
    selected: Optional[List[str]] = (
        [f.strip() for f in fields.split(",")] if fields else None
    )
    logger.info("Building Ensembl gene table")
    _create_ensembl_gene_tb(
        input_path=raw_input,
        output_path=output_ht,
        fields=selected,
        canonical=canonical,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"Ensembl gene table created at {output_ht}")


@mktable_group.command("dbnsfp")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
@_ref_genome_opt
@click.option(
    "--no-parse-transcript-scores",
    is_flag=True,
    help="Do not parse transcript-specific score fields into dicts",
)
@click.option(
    "--group-prefixes",
    type=str,
    default=None,
    help="Comma-separated prefixes to group into structs (default: gnomAD,ExAC,1000Gp3,ESP6500,clinvar)",
)
def mktable_dbnsfp(
    raw_input: str,
    output_ht: str,
    overwrite: bool,
    export_tsv: bool,
    ref_genome: str,
    no_parse_transcript_scores: bool,
    group_prefixes: Optional[str],
):
    """Build a dbNSFP variant annotation Table from a TSV/BGZ (keyed by locus, alleles)."""
    prefixes = None
    if group_prefixes:
        prefixes = [p.strip() for p in group_prefixes.split(",") if p.strip()]

    logger.info("Building dbNSFP table")
    _create_dbnsfp_tb(
        input_path=raw_input,
        output_path=output_ht,
        reference_genome=ref_genome,
        overwrite=overwrite,
        export_tsv=export_tsv,
        parse_transcript_scores=not no_parse_transcript_scores,
        group_prefixes=prefixes,
    )
    click.echo(f"dbNSFP table created at {output_ht}")


@mktable_group.command("clingen-gene-disease")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
@click.option(
    "--key-by",
    type=click.Choice(["gene_disease", "gene"], case_sensitive=False),
    default="gene_disease",
    show_default=True,
    help="Keying strategy: 'gene_disease' (hgnc_id, mondo_id) or 'gene' (aggregated by hgnc_id)",
)
@click.option(
    "--min-classification",
    type=click.Choice(
        ["Definitive", "Strong", "Moderate", "Limited", "Disputed", "Refuted"],
        case_sensitive=True,
    ),
    default=None,
    help="Filter to classifications at or above this level",
)
@click.option(
    "--fields",
    type=str,
    default=None,
    help="Comma-separated list of fields to retain (optional)",
)
def mktable_clingen_gene_disease(
    raw_input: str,
    output_ht: str,
    overwrite: bool,
    export_tsv: bool,
    key_by: str,
    min_classification: Optional[str],
    fields: Optional[str],
):
    """Build a ClinGen Gene-Disease Validity Table from a CSV (keyed by gene or gene-disease)."""
    selected: Optional[List[str]] = (
        [f.strip() for f in fields.split(",")] if fields else None
    )
    logger.info("Building ClinGen Gene-Disease table")
    _create_clingen_gene_disease_tb(
        input_path=raw_input,
        output_path=output_ht,
        key_by=key_by.lower(),
        min_classification=min_classification,
        fields=selected,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"ClinGen Gene-Disease table created at {output_ht}")


@mktable_group.command("hgnc")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
@click.option(
    "--include-withdrawn",
    is_flag=True,
    help="Include withdrawn/non-approved genes (default: only approved)",
)
@click.option(
    "--fields",
    type=str,
    default=None,
    help="Comma-separated list of fields to retain (optional)",
)
def mktable_hgnc(
    raw_input: str,
    output_ht: str,
    overwrite: bool,
    export_tsv: bool,
    include_withdrawn: bool,
    fields: Optional[str],
):
    """Build an HGNC gene nomenclature Table from a TSV (keyed by hgnc_id)."""
    selected: Optional[List[str]] = (
        [f.strip() for f in fields.split(",")] if fields else None
    )
    logger.info("Building HGNC gene table")
    _create_hgnc_gene_tb(
        input_path=raw_input,
        output_path=output_ht,
        include_withdrawn=include_withdrawn,
        fields=selected,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"HGNC gene table created at {output_ht}")
