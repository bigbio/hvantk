"""
A module to create Hail Tables from different raw data files.
The tables are used to annotate variants and genes with different features.

"""

import hail as hl
import logging

logger = logging.getLogger(__name__)

from hvantk.utils.constants import ENSEMBL_BIOMART_FIELDS

from hvantk.settings import RAW_DATA_PATHS

raw_resource_paths = RAW_DATA_PATHS


def create_gnomad_constraint_gene_metrics_tb(
    input_path: str,
    output_path: str,
    fields: list = None,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> hl.Table:
    """
    Creates a Hail Table of gene-level constraint metrics from a gnomAD input file.
    
    Imports the table keyed by gene ID, optionally selects specified fields, checkpoints
    the result to disk, and can export it as a compressed TSV file.
    
    Args:
        input_path: Path to the gnomAD constraint metrics file.
        output_path: Destination path for the checkpointed Hail Table.
        fields: Optional list of fields to retain in the table.
        overwrite: If True, overwrites any existing output at the destination.
        export_tsv: If True, exports the table as a compressed TSV file.
    
    Returns:
        A Hail Table containing gene-level constraint metrics.
    """
    logger.info(f"Creating gnomAD constraint gene metrics table from {input_path}")
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
) -> hl.Table:
    """
    Creates a Hail Table of protein-protein interactions from a BED file.
    
    Imports interaction intervals, removes duplicates, and checkpoints the resulting table to disk. Optionally exports the table as a compressed TSV file.
    
    Args:
        input_path: Path to the input BED file containing interaction data.
        output_path: Destination path for the checkpointed Hail Table.
        overwrite: If True, overwrites any existing output at the destination.
        export_tsv: If True, exports the table as a compressed TSV file.
        reference_genome: Reference genome build to use for interval parsing.
    
    Returns:
        A Hail Table containing the processed protein-protein interaction data.
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
) -> hl.Table:
    """
    Creates a Hail Table of ClinVar variants from a VCF file.
    
    Imports ClinVar variant data, recodes contig names to the "chr" format, skips invalid loci, and keys the table by locus and alleles. The resulting table is checkpointed to disk and can optionally be exported as a compressed TSV file.
    
    Args:
        input_path: Path to the input ClinVar VCF file.
        output_path: Destination path for the checkpointed Hail Table.
        overwrite: If True, overwrites any existing output at the destination.
        export_tsv: If True, exports the table as a compressed TSV file.
        reference_genome: Reference genome to use for import (default: "GRCh38").
    
    Returns:
        A Hail Table containing ClinVar variant annotations.
    """
    logger.info(f"Creating ClinVar table from {input_path}")
    recode = {f"{i}": f"chr{i}" for i in (list(range(1, 23)) + ["X", "Y"])}
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
) -> hl.Table:
    """
    Creates a Hail Table of gene-level constraint metrics from a GeVir input file.
    
    Args:
        input_path: Path to the GeVir input file.
        output_path: Destination path for the checkpointed Hail Table.
        fields: Optional list of fields to select from the imported table.
        overwrite: If True, overwrites any existing output at the destination.
        export_tsv: If True, exports the resulting table as a compressed TSV file.
    
    Returns:
        A Hail Table keyed by gene ID containing GeVir constraint metrics.
    """
    logger.info(f"Creating GEVIR table from {input_path}")
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
) -> hl.Table:
    """
    Creates a Hail Table of gene-level metrics from an Ensembl Biomart export.
    
    Aggregates transcript, protein, and synonym information per gene, with optional filtering for canonical transcripts and selection of specific fields. The resulting table is checkpointed to disk and can be exported as a TSV file.
    
    Args:
        input_path: Path to the Ensembl Biomart export file.
        output_path: Destination path for the checkpointed Hail Table.
        fields: Optional list of fields to include in the output table.
        canonical: If True, includes only canonical transcripts.
        overwrite: If True, overwrites any existing output file.
        export_tsv: If True, exports the table as a compressed TSV file.
    
    Returns:
        A Hail Table keyed by gene ID with aggregated gene-level metrics.
    """
    logger.info(f"Creating Ensembl gene table from {input_path}")
    gene_tb = hl.import_table(paths=input_path, min_partitions=50)

    # replace field names
    logger.info("Replacing field names")
    gene_tb = gene_tb.rename(ENSEMBL_BIOMART_FIELDS)

    # filter canonical transcripts
    if canonical:
        logger.info("Filtering canonical transcripts")
        gene_tb = gene_tb.filter(gene_tb.canonical == "1", keep=True)

    # group by gene_id
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
            gene_name=hl.agg.collect(gene_tb.gene_name).first(),
            chromosome=hl.agg.collect(gene_tb.chromosome).first(),
            gene_start=hl.agg.collect(gene_tb.gene_start).first(),
            gene_end=hl.agg.collect(gene_tb.gene_end).first(),
            gene_type=hl.agg.collect(gene_tb.gene_type).first(),
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
