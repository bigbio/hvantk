"""
Hail Table builders for converting raw sources into Hail Tables (HT).

This module supersedes 'creators.py'. Prefer importing from 'table_builders'.
"""

import hail as hl
import logging
from typing import Optional, List
from hvantk.utils.table_utils import get_row_fields

logger = logging.getLogger(__name__)

from hvantk.core.constants import ENSEMBL_BIOMART_FIELDS
from hvantk.utils.genome import contig_recoding  # correct module import

__all__ = [
    "create_gnomad_constraint_gene_metrics_tb",
    "create_interactome_tb",
    "create_clinvar_tb",
    "create_gevir_tb",
    "create_ensembl_gene_tb",
    "create_dbnsfp_tb",
]


def create_gnomad_constraint_gene_metrics_tb(
    input_path: str,
    output_path: str,
    fields: list = None,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> 'hl.Table':
    """
    Create a Hail Table from gnomAD constraint gene metrics TSV file keyed by gene_id.
    Example usage:
        ht = create_gnomad_constraint_gene_metrics_tb(
            input_path="/path/to/gnomad_constraint_metrics.txt",
            output_path="/path/to/output.ht",
            fields=["gene_id", "pLI", "oe_lof"]
        )
    Parameters
    ----------
    input_path : str
        Path to the gnomAD constraint gene metrics TSV input file.
    output_path : str
        Path to write the output Hail Table.
    fields : list, optional
        List of fields to select from the imported table (default: None, selects all fields).
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).
    Returns
    -------
    hl.Table
        Hail Table keyed by gene_id with selected constraint metrics.
    """

    logger.info(f"Creating gnomAD constraint gene metrics table from {input_path}")
    # If schema is stable, consider specifying types=... instead of impute=True
    gnomad_tb = hl.import_table(
        paths=input_path, impute=True, min_partitions=100, key="gene_id"
    )

    if fields is not None:
        logger.info(f"Selecting fields: {fields}")
        gnomad_tb = gnomad_tb.select(*fields)

    logger.info(f"Checkpointing table to {output_path}")
    gnomad_tb = gnomad_tb.checkpoint(output=output_path, overwrite=overwrite)

    if export_tsv:
        logger.info(f"Exporting table to {output_path}.tsv.bgz")
        gnomad_tb.export(output_path + ".tsv.bgz")

    return gnomad_tb


def create_interactome_tb(
    input_path: str,
    output_path: str,
    overwrite: bool = False,
    export_tsv: bool = False,
    reference_genome: str = "GRCh38",
) -> 'hl.Table':
    """
    Create a Hail Table from a protein-protein interaction BED file.
    Example usage:
        ht = create_interactome_tb(
            input_path="/path/to/interactome.bed",
            output_path="/path/to/output.ht"
        )
    Parameters
    ----------
    input_path : str
        Path to the protein-protein interaction BED input file.
    output_path : str
        Path to write the output Hail Table.
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).
    reference_genome : str, optional
        Reference genome to use for parsing intervals (default: "GRCh38").
    Returns
    -------
    hl.Table
        Hail Table with protein-protein interactions.
    """

    logger.info(f"Creating interactome table from {input_path}")
    ppi_tb = (
        hl.import_bed(
            path=input_path,
            skip_invalid_intervals=True,
            reference_genome=reference_genome,
        )
        .repartition(100)
        .distinct()
    )

    logger.info(f"Checkpointing table to {output_path}")
    ppi_tb = ppi_tb.checkpoint(output=output_path, overwrite=overwrite)

    if export_tsv:
        logger.info(f"Exporting table to {output_path}.tsv.bgz")
        ppi_tb.export(f"{output_path}.tsv.bgz")

    return ppi_tb


def create_clinvar_tb(
    input_path: str,
    output_path: str,
    overwrite: bool = False,
    export_tsv: bool = False,
    reference_genome: str = "GRCh38",
) -> 'hl.Table':
    """
    Create a Hail Table from a ClinVar VCF file keyed by (locus, alleles).
    Example usage:
        ht = create_clinvar_tb(
            input_path="/path/to/clinvar.vcf.gz",
            output_path="/path/to/output.ht"
        )
    Parameters
    ----------
    input_path : str
        Path to the ClinVar VCF input file.
    output_path : str
        Path to write the output Hail Table.
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a flattened TSV version (default: False).
    reference_genome : str, optional
        Reference genome to use for parsing variants (default: "GRCh38").
    Returns
    -------
    hl.Table
        Hail Table keyed by (locus, alleles) with ClinVar annotations.
    """
    logger.info(f"Creating ClinVar table from {input_path}")
    # Use utility for contig recoding
    recode = contig_recoding()
    clinvar_tb = (
        hl.import_vcf(
            path=input_path,
            force=True,
            reference_genome=reference_genome,
            contig_recoding=recode,
            skip_invalid_loci=True,
        )
        .rows()
        .repartition(100)
        .key_by("locus", "alleles")
    )

    logger.info(f"Checkpointing table to {output_path}")
    clinvar_tb = clinvar_tb.checkpoint(output=output_path, overwrite=overwrite)

    if export_tsv:
        logger.info(f"Exporting table to {output_path}.tsv.bgz")
        (clinvar_tb.flatten().export(f"{output_path}.tsv.bgz"))

    return clinvar_tb


def create_gevir_tb(
    input_path: str,
    output_path: str,
    fields: list = None,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> 'hl.Table':
    """
    Create a Hail Table from GEVIR gene metrics TSV file keyed by gene_id.
    Example usage:
        ht = create_gevir_tb(
            input_path="/path/to/gevir_metrics.txt",
            output_path="/path/to/output.ht",
            fields=["gene_id", "gevir_score", "gevir_rank"]
        )
    Parameters
    ----------
    input_path : str
        Path to the GEVIR gene metrics TSV input file.
    output_path : str
        Path to write the output Hail Table.
    fields : list, optional
        List of fields to select from the imported table (default: None, selects all fields).
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).
    Returns
    -------
    hl.Table
        Hail Table keyed by gene_id with selected GEVIR metrics.

    """
    logger.info(f"Creating GEVIR table from {input_path}")
    # If schema is stable, consider specifying types=... instead of impute=True
    gevir_tb = hl.import_table(
        paths=input_path, impute=True, min_partitions=100, key="gene_id"
    )

    if fields is not None:
        logger.info(f"Selecting fields: {fields}")
        gevir_tb = gevir_tb.select(*fields)

    logger.info(f"Checkpointing table to {output_path}")
    gevir_tb = gevir_tb.checkpoint(output=output_path, overwrite=overwrite)

    if export_tsv:
        logger.info(f"Exporting table to {output_path}.tsv.bgz")
        gevir_tb.export(output_path + ".tsv.bgz")

    return gevir_tb


def create_ensembl_gene_tb(
    input_path: str,
    output_path: str,
    fields: list = None,
    canonical: bool = True,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> 'hl.Table':
    """
    Create a Hail Table from an Ensembl BioMart gene annotation TSV file keyed by gene_id.
    Example usage:
        ht = create_ensembl_gene_tb(
            input_path="/path/to/ensembl_biomart_genes.txt",
            output_path="/path/to/output.ht",
            fields=["gene_id", "gene_name", "gene_type", "chromosome",
                    "gene_start", "gene_end", "transcript_id", "protein_id"]
        )
    Parameters
    ----------
    input_path : str
        Path to the Ensembl BioMart gene annotation TSV input file.
    output_path : str
        Path to write the output Hail Table.
    fields : list, optional
        List of fields to select from the imported table (default: None, selects all fields).
    canonical : bool, optional
        If True, filter to only canonical transcripts (default: True).
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).
    Returns
    -------
    hl.Table
        Hail Table keyed by gene_id with selected Ensembl gene annotations.
    """
    logger.info(f"Creating Ensembl gene table from {input_path}")
    gene_tb = hl.import_table(paths=input_path, min_partitions=50, impute=True)

    logger.info("Replacing field names")
    gene_tb = gene_tb.rename(ENSEMBL_BIOMART_FIELDS)

    if canonical:
        logger.info("Filtering canonical transcripts")
        gene_tb = gene_tb.filter(gene_tb.canonical == "1", keep=True)

    logger.info("Grouping by gene_id")
    gene_tb = (
        gene_tb.group_by(gene_tb.gene_id)
        .aggregate(
            transcript_id=hl.agg.collect_as_set(gene_tb.transcript_id).filter(
                lambda x: x != ""
            ),
            protein_id=hl.agg.collect_as_set(gene_tb.protein_id).filter(
                lambda x: x != ""
            ),
            gene_synonym=hl.agg.collect_as_set(gene_tb.gene_synonym).filter(
                lambda x: x != ""
            ),
            gene_name=hl.agg.take(hl.or_else(gene_tb.gene_name, ""), 1)[0],
            chromosome=hl.agg.take(hl.or_else(gene_tb.chromosome, ""), 1)[0],
            gene_start=hl.agg.take(hl.or_else(gene_tb.gene_start, hl.null(gene_tb.gene_start.dtype)), 1)[0],
            gene_end=hl.agg.take(hl.or_else(gene_tb.gene_end, hl.null(gene_tb.gene_end.dtype)), 1)[0],
            gene_type=hl.agg.take(hl.or_else(gene_tb.gene_type, ""), 1)[0],
        )
        .key_by("gene_id")
    )

    if fields is not None:
        logger.info(f"Selecting fields: {fields}")
        gene_tb = gene_tb.select(*fields)

    logger.info(f"Checkpointing table to {output_path}")
    gene_tb = gene_tb.checkpoint(output=output_path, overwrite=overwrite)

    if export_tsv:
        logger.info(f"Exporting table to {output_path}.tsv.bgz")
        gene_tb.export(output_path + ".tsv.bgz")

    return gene_tb


def create_dbnsfp_tb(
    input_path: str,
    output_path: str,
    reference_genome: str = "GRCh38",
    overwrite: bool = False,
    export_tsv: bool = False,
    min_partitions: int = 200,
    force_bgz: bool = True,
    parse_transcript_scores: bool = True,
    group_prefixes: Optional[List[str]] = None,
) -> 'hl.Table':
    """
    Create a Hail Table from a dbNSFP variant TSV/BGZ file keyed by (locus, alleles).

    Example usage:
        ht = create_dbnsfp_tb(
            input_path="/path/to/dbNSFP.tsv.bgz",
            output_path="/path/to/output.ht"
        )

    Parameters
    ----------
    input_path : str
        Path to the dbNSFP TSV/BGZ input file.
    output_path : str
        Path to write the output Hail Table.
    reference_genome : str, optional
        Reference genome to use for parsing variants (default: "GRCh38").
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a flattened TSV version (default: False).
    min_partitions : int, optional
        Minimum number of partitions for import (default: 200).
    force_bgz : bool, optional
        If True, force bgzip compression for input (default: True).
    parse_transcript_scores : bool, optional
        If True, parse transcript-specific scores into dicts (default: True).
    group_prefixes : list of str, optional
        List of field prefixes to group into structs (default: common population/annotation prefixes).

    Returns
    -------
    hl.Table
        Hail Table keyed by (locus, alleles) with parsed and grouped annotations.

    Steps performed:
    - Import table with missing '.' and no type imputation
    - Build a variant key from '#chr', 'pos(1-based)', 'ref', 'alt' and parse to (locus, alleles)
    - Key the table by (locus, alleles)
    - Optionally map transcript-specific scores ending with '_score' or 'CADD_phred' to dict(Ensembl_transcriptid -> float)
    - Optionally group common prefixes (e.g., gnomAD, ExAC) into structs and drop original prefixed columns
    """

    logger.info(f"Importing dbNSFP table from {input_path}")
    ht = hl.import_table(
        paths=input_path,
        min_partitions=min_partitions,
        impute=False,
        missing='.',
        force_bgz=force_bgz,
    )

    # Normalize chromosome field and construct variant key
    row_fields = get_row_fields(ht)
    if "#chr" in row_fields:
        ht = ht.rename({'#chr': 'chr'})
    else:
        # Some exports might already use 'chr'
        if "chr" not in row_fields:
            raise ValueError("dbNSFP input missing '#chr' or 'chr' column")

    _chr_str = hl.str(ht['chr'])
    ht = ht.annotate(
        chr=hl.if_else(_chr_str.lower().startswith("chr"), _chr_str, hl.str("chr") + _chr_str)
    )

    # Build variant_key: chr:pos:ref:alt
    row_fields = get_row_fields(ht)
    if 'pos(1-based)' not in row_fields or 'ref' not in row_fields or 'alt' not in row_fields:
        raise ValueError("dbNSFP input missing required columns: 'pos(1-based)', 'ref', or 'alt'")

    variant_key_expr = hl.array([
        ht.chr,
        hl.str(ht['pos(1-based)']),
        ht.ref,
        ht.alt,
    ])
    ht = ht.annotate(variant_key=hl.delimit(variant_key_expr, ':'))

    # Parse to locus/alleles
    ht = ht.annotate(**hl.parse_variant(ht.variant_key, reference_genome=reference_genome))

    # Key the table by (locus, alleles) before any selects to avoid overwriting key fields
    ht = ht.key_by('locus', 'alleles')
    # Optional cleanup of staging columns; keep if downstream needs them
    ht = ht.drop('variant_key')
    ht = ht.drop('chr', 'pos(1-based)', 'ref', 'alt')

    # Transcript-specific score parsing
    row_fields = get_row_fields(ht)
    if parse_transcript_scores and 'Ensembl_transcriptid' in row_fields:
        logger.info("Parsing transcript-specific scores into dicts keyed by Ensembl_transcriptid")
        ht = ht.annotate(Ensembl_transcriptid=hl.str(ht.Ensembl_transcriptid))
        ht = ht.annotate(Ensembl_transcriptid=ht.Ensembl_transcriptid.split(";"))

        row_fields_list = list(get_row_fields(ht))
        score_fields = [f for f in row_fields_list if f.endswith('_score') or f == 'CADD_phred']
        def _to_float_array(s):
            s_def = hl.or_else(s, "")  # empty string if missing
            arr = s_def.split(";")
            return hl.map(lambda x: hl.parse_float(x), arr)

        def _single_to_dict(val):
            # Map same scalar value to all transcripts
            return hl.dict(hl.zip(ht.Ensembl_transcriptid,
                                  hl.map(lambda _x: hl.parse_float(val), ht.Ensembl_transcriptid)))

        ann = {}
        for f in score_fields:
            is_multi = hl.is_defined(ht[f]) & ht[f].contains(";")
            ann[f] = hl.if_else(
                is_multi,
                hl.dict(hl.zip(ht.Ensembl_transcriptid, _to_float_array(ht[f]))),
                _single_to_dict(ht[f])
            )
        if ann:
            ht = ht.annotate(**ann)

    # Group common prefixes into structs and drop original columns
    if group_prefixes is None:
        group_prefixes = ['gnomAD', 'ExAC', '1000Gp3', 'ESP6500', 'clinvar']

    for prefix in group_prefixes:
        row_fields_list = list(get_row_fields(ht))
        pref_fields = [f for f in row_fields_list if f != prefix and f.startswith(prefix)]
        if pref_fields:
            logger.info(f"Grouping {prefix}* fields into struct '{prefix}' ({len(pref_fields)} fields)")
            ht = ht.annotate(**{prefix: hl.struct(**{f: ht[f] for f in pref_fields})})
            # Drop original columns (preserve keys implicitly)
            ht = ht.drop(*pref_fields)

    logger.info(f"Checkpointing dbNSFP table to {output_path}")
    ht = ht.checkpoint(output=output_path, overwrite=overwrite)

    if export_tsv:
        logger.info(f"Exporting flattened dbNSFP table to {output_path}.tsv.bgz")
        ht.flatten().export(f"{output_path}.tsv.bgz")

    return ht
