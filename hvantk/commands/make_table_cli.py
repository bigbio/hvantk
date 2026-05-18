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
    from hvantk.skills.clinvar.builder import create_clinvar_tb

    return create_clinvar_tb(*args, **kwargs)


def _create_interactome_tb(*args, **kwargs):
    from hvantk.skills.insider.builder import create_interactome_tb

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
    from hvantk.skills.clingen.builder import create_clingen_gene_disease_tb

    return create_clingen_gene_disease_tb(*args, **kwargs)


def _create_gencc_submissions_tb(*args, **kwargs):
    from hvantk.skills.gencc.builder import create_gencc_submissions_tb

    return create_gencc_submissions_tb(*args, **kwargs)


def _create_hgnc_gene_tb(*args, **kwargs):
    from hvantk.skills.hgnc.builder import create_hgnc_gene_tb

    return create_hgnc_gene_tb(*args, **kwargs)


def _create_ptm_sites_tb(*args, **kwargs):
    from hvantk.skills.uniprot_ptm.builder import create_ptm_sites_tb

    return create_ptm_sites_tb(*args, **kwargs)


def _create_eqtl_tb(*args, **kwargs):
    from hvantk.skills.gtex_eqtl.builder import create_eqtl_tb

    return create_eqtl_tb(*args, **kwargs)


def _create_pqtl_tb(*args, **kwargs):
    from hvantk.tables.table_builders import create_pqtl_tb

    return create_pqtl_tb(*args, **kwargs)


def _create_gwas_catalog_tb(*args, **kwargs):
    from hvantk.skills.gwas_catalog.builder import create_gwas_catalog_tb

    return create_gwas_catalog_tb(*args, **kwargs)


def _create_msigdb_tb(*args, **kwargs):
    from hvantk.skills.msigdb.builder import create_msigdb_tb

    return create_msigdb_tb(*args, **kwargs)


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
    ht = _create_clinvar_tb(
        input_path=raw_input,
        output_path=output_ht,
        overwrite=overwrite,
        export_tsv=export_tsv,
        reference_genome=ref_genome,
    )
    click.echo(f"ClinVar table created at {output_ht}")
    ht.describe()


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
    ht = _create_interactome_tb(
        input_path=raw_input,
        output_path=output_ht,
        overwrite=overwrite,
        export_tsv=export_tsv,
        reference_genome=ref_genome,
    )
    click.echo(f"Interactome table created at {output_ht}")
    ht.describe()


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
    ht = _create_gevir_tb(
        input_path=raw_input,
        output_path=output_ht,
        fields=selected,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"GEVIR table created at {output_ht}")
    ht.describe()


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
    ht = _create_gnomad_constraint_gene_metrics_tb(
        input_path=raw_input,
        output_path=output_ht,
        fields=selected,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"gnomAD metrics table created at {output_ht}")
    ht.describe()


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
    ht = _create_ensembl_gene_tb(
        input_path=raw_input,
        output_path=output_ht,
        fields=selected,
        canonical=canonical,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"Ensembl gene table created at {output_ht}")
    ht.describe()


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
@click.option(
    "--auto-convert-bgz",
    is_flag=True,
    help="Automatically convert plain gzip files to BGZF before import",
)
def mktable_dbnsfp(
    raw_input: str,
    output_ht: str,
    overwrite: bool,
    export_tsv: bool,
    ref_genome: str,
    no_parse_transcript_scores: bool,
    group_prefixes: Optional[str],
    auto_convert_bgz: bool,
):
    """Build a dbNSFP variant annotation Table from a TSV/BGZ (keyed by locus, alleles)."""
    prefixes = None
    if group_prefixes:
        prefixes = [p.strip() for p in group_prefixes.split(",") if p.strip()]

    logger.info("Building dbNSFP table")
    ht = _create_dbnsfp_tb(
        input_path=raw_input,
        output_path=output_ht,
        reference_genome=ref_genome,
        overwrite=overwrite,
        export_tsv=export_tsv,
        parse_transcript_scores=not no_parse_transcript_scores,
        group_prefixes=prefixes,
        auto_convert_bgz=auto_convert_bgz,
    )
    click.echo(f"dbNSFP table created at {output_ht}")
    ht.describe()


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
    ht = _create_clingen_gene_disease_tb(
        input_path=raw_input,
        output_path=output_ht,
        key_by=key_by.lower(),
        min_classification=min_classification,
        fields=selected,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"ClinGen Gene-Disease table created at {output_ht}")
    ht.describe()


@mktable_group.command("gencc-submissions")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
@click.option(
    "--key-by",
    type=click.Choice(
        ["gene_disease_submitter", "gene_disease", "gene"], case_sensitive=False
    ),
    default="gene_disease_submitter",
    show_default=True,
    help=(
        "Keying strategy: 'gene_disease_submitter' (full granularity), "
        "'gene_disease' (aggregate submitters), or 'gene' (aggregate by gene)"
    ),
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
        ],
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
def mktable_gencc_submissions(
    raw_input: str,
    output_ht: str,
    overwrite: bool,
    export_tsv: bool,
    key_by: str,
    min_classification: Optional[str],
    fields: Optional[str],
):
    """Build a GenCC Submissions Table from a TSV (keyed by gene, gene-disease, or gene-disease-submitter)."""
    selected: Optional[List[str]] = (
        [f.strip() for f in fields.split(",")] if fields else None
    )
    logger.info("Building GenCC Submissions table")
    ht = _create_gencc_submissions_tb(
        input_path=raw_input,
        output_path=output_ht,
        key_by=key_by.lower(),
        min_classification=min_classification,
        fields=selected,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"GenCC Submissions table created at {output_ht}")
    ht.describe()


@mktable_group.command("cosmic-cgc")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
@click.option(
    "--hgnc-path",
    type=str,
    default=None,
    help="Path to HGNC Hail Table (.ht) for gene symbol -> HGNC ID resolution",
)
@click.option(
    "--mutation-context",
    type=click.Choice(["somatic", "germline", "both"]),
    default="both",
    show_default=True,
    help="Filter genes by somatic/germline mutation context",
)
@click.option(
    "--min-classification",
    type=click.Choice(["Tier 1", "Tier 2"], case_sensitive=True),
    default=None,
    help="Filter to classifications at or above this level",
)
def mktable_cosmic_cgc(
    raw_input: str,
    output_ht: str,
    overwrite: bool,
    export_tsv: bool,
    hgnc_path: Optional[str],
    mutation_context: str,
    min_classification: Optional[str],
):
    """Build a COSMIC Cancer Gene Census Hail Table from the downloaded TSV."""
    logger.info("Building COSMIC CGC table")
    from hvantk.tables.table_builders import create_cosmic_cgc_tb

    ht = create_cosmic_cgc_tb(
        input_path=raw_input,
        output_path=output_ht,
        hgnc_path=hgnc_path,
        min_classification=min_classification,
        mutation_context=mutation_context,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"COSMIC CGC table created at {output_ht}")
    ht.describe()


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
    ht = _create_hgnc_gene_tb(
        input_path=raw_input,
        output_path=output_ht,
        include_withdrawn=include_withdrawn,
        fields=selected,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"HGNC gene table created at {output_ht}")
    ht.describe()


@mktable_group.command("ptm-sites")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
@_ref_genome_opt
@click.option(
    "--flanking-codons",
    type=int,
    default=5,
    show_default=True,
    help="Number of flanking codons for proximal window",
)
def mktable_ptm_sites(
    raw_input: str,
    output_ht: str,
    overwrite: bool,
    export_tsv: bool,
    ref_genome: str,
    flanking_codons: int,
):
    """Build a PTM sites Hail Table from mapped coordinates (keyed by locus)."""
    logger.info("Building PTM sites table")
    ht = _create_ptm_sites_tb(
        input_path=raw_input,
        output_path=output_ht,
        reference_genome=ref_genome,
        flanking_codons=flanking_codons,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"PTM sites table created at {output_ht}")
    ht.describe()


# ---------------------------------------------------------------------------
# eQTL builder
# ---------------------------------------------------------------------------


@mktable_group.command("eqtl")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
@_ref_genome_opt
@click.option(
    "--source",
    type=click.Choice(["gtex_v11", "gtex_v8", "eqtlgen"]),
    default="gtex_v11",
    show_default=True,
    help="eQTL data source format",
)
@click.option(
    "--tissue",
    type=str,
    default=None,
    help="Restrict import to files matching this tissue name",
)
@click.option(
    "--p-threshold",
    type=float,
    default=5e-8,
    show_default=True,
    help="P-value threshold (0 to keep all pairs for coloc)",
)
def mktable_eqtl(
    raw_input: str,
    output_ht: str,
    overwrite: bool,
    export_tsv: bool,
    ref_genome: str,
    source: str,
    tissue: str,
    p_threshold: float,
):
    """Build an eQTL Hail Table (keyed by locus, alleles, gene_id)."""
    logger.info("Building eQTL table (source=%s)", source)
    ht = _create_eqtl_tb(
        input_path=raw_input,
        output_path=output_ht,
        reference_genome=ref_genome,
        source=source,
        tissue=tissue,
        p_threshold=p_threshold,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"eQTL table created at {output_ht}")
    ht.describe()


# ---------------------------------------------------------------------------
# pQTL builder
# ---------------------------------------------------------------------------


@mktable_group.command("pqtl")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
@_ref_genome_opt
@click.option(
    "--source",
    type=click.Choice(["gtex_fang"]),
    default="gtex_fang",
    show_default=True,
    help="pQTL data source format",
)
@click.option("--tissue", type=str, default=None, help="Restrict to this tissue")
@click.option(
    "--hgnc-ht",
    type=str,
    default=None,
    help="HGNC Hail Table (built by 'hvantk mktable hgnc-gene') for "
    "gene symbol → Ensembl ID mapping. Required unless --no-gene-map.",
)
@click.option(
    "--no-gene-map",
    is_flag=True,
    default=False,
    help="Skip Ensembl mapping — key by raw gene symbol. "
    "The table will NOT join with eQTL tables in cascade analysis.",
)
@click.option(
    "--p-threshold",
    type=float,
    default=None,
    help="P-value threshold (omit to keep all pairs)",
)
def mktable_pqtl(
    raw_input: str,
    output_ht: str,
    overwrite: bool,
    export_tsv: bool,
    ref_genome: str,
    source: str,
    tissue: str,
    hgnc_ht: str,
    no_gene_map: bool,
    p_threshold: float,
):
    """Build a pQTL Hail Table (keyed by locus, alleles, gene_id)."""
    logger.info("Building pQTL table (source=%s)", source)
    ht = _create_pqtl_tb(
        input_path=raw_input,
        output_path=output_ht,
        reference_genome=ref_genome,
        source=source,
        tissue=tissue,
        hgnc_ht=hgnc_ht,
        no_gene_map=no_gene_map,
        p_threshold=p_threshold,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"pQTL table created at {output_ht}")
    ht.describe()


# ---------------------------------------------------------------------------
# AlphaGenome builder
# ---------------------------------------------------------------------------


@mktable_group.command("alphagenome")
@click.option(
    "--input",
    "input_path",
    required=True,
    type=str,
    help="Path to Hail Table (.ht) or TSV with chrom/pos/ref/alt columns",
)
@click.option(
    "--output-dir",
    required=True,
    type=str,
    help="Output directory for per-modality Hail Tables",
)
@click.option(
    "--config",
    "config_path",
    required=True,
    type=str,
    help="Path to AlphaGenome YAML config file",
)
@click.option(
    "--no-resume",
    is_flag=True,
    help="Discard existing checkpoints and restart from scratch",
)
@_overwrite_opt
def mktable_alphagenome(input_path, output_dir, config_path, no_resume, overwrite):
    """Run AlphaGenome variant effect predictions."""
    from hvantk.tables.table_builders import create_alphagenome_tb

    logger.info("Running AlphaGenome variant predictions")
    create_alphagenome_tb(
        input_path=input_path,
        output_path=output_dir,
        config_path=config_path,
        no_resume=no_resume,
        overwrite=overwrite,
    )
    click.echo(f"AlphaGenome predictions written to {output_dir}")


@mktable_group.command("gwas-catalog")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
@_ref_genome_opt
def mktable_gwas_catalog(
    raw_input: str, output_ht: str, overwrite: bool, export_tsv: bool, ref_genome: str
):
    """Build a GWAS Catalog Hail Table from the v1.0 full-associations TSV (keyed by locus, alleles)."""
    logger.info("Building GWAS Catalog table")
    ht = _create_gwas_catalog_tb(
        input_path=raw_input,
        output_path=output_ht,
        overwrite=overwrite,
        export_tsv=export_tsv,
        reference_genome=ref_genome,
    )
    click.echo(f"GWAS Catalog table created at {output_ht}")
    ht.describe()


@mktable_group.command("msigdb")
@_raw_input_opt
@_output_ht_opt
@_overwrite_opt
@_export_tsv_opt
def mktable_msigdb(
    raw_input: str, output_ht: str, overwrite: bool, export_tsv: bool
):
    """Build an MSigDB Hail Table from a GMT file (keyed by set_name)."""
    logger.info("Building MSigDB table")
    ht = _create_msigdb_tb(
        input_path=raw_input,
        output_path=output_ht,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
    click.echo(f"MSigDB table created at {output_ht}")
    ht.describe()
