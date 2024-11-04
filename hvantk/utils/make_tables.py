"""
A module to create Hail Tables from different raw data files.
The tables are used to annotate variants and genes with different features.

"""

import hail as hl

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
    Create a Hail Table with gene-level constraint metrics from gnomad.

    :param input_path: Path to the input file.
    :param output_path: Path to save the output Hail Table.
    :param fields: List of fields to select from the input table.
    :param overwrite: Whether to overwrite the existing output file.
    :param export_tsv: Whether to export the table as a TSV file.
    :return: Hail Table
    """
    gnomad_tb = hl.import_table(paths=input_path,
                                impute=True,
                                min_partitions=100,
                                key='gene_id')

    if fields is not None:
        gnomad_tb = gnomad_tb.select(*fields)

    gnomad_tb = gnomad_tb.checkpoint(output=output_path, overwrite=overwrite)

    if export_tsv:
        gnomad_tb.export(output_path + '.tsv.bgz')

    return gnomad_tb


def create_interactome_tb(
        input_path: str,
        output_path: str,
        overwrite: bool = False,
        export_tsv: bool = False,
        reference_genome: str = 'GRCh38',
) -> hl.Table:
    """
    Create a Hail Table with protein-protein interaction data.

    :param input_path: Path to the input BED file.
    :param output_path: Path to save the output Hail Table.
    :param overwrite: Whether to overwrite the existing output file.
    :param export_tsv: Whether to export the table as a TSV file.
    :param reference_genome: Reference genome to use (default: GRCh38).

    :return: Hail Table
    """
    ppi_tb = hl.import_bed(
        path=input_path,
        skip_invalid_intervals=True,
        reference_genome=reference_genome
    ).repartition(100).distinct()

    ppi_tb = ppi_tb.checkpoint(output=output_path, overwrite=overwrite)

    if export_tsv:
        ppi_tb.export(f"{output_path}.tsv.bgz")

    return ppi_tb


def create_clinvar_tb(
        input_path: str,
        output_path: str,
        overwrite: bool = False,
        export_tsv: bool = False,
        reference_genome: str = 'GRCh38'
) -> hl.Table:
    """
    Create a Hail Table with clinvar variants.

    :param input_path: Path to the input VCF file.
    :param output_path: Path to save the output Hail Table.
    :param overwrite: Whether to overwrite the existing output file.
    :param export_tsv: Whether to export the table as a TSV file.
    :param reference_genome: Reference genome to use (default: GRCh38).
    :return: Hail Table
    """
    recode = {f"{i}": f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
    clinvar_tb = (hl.import_vcf(path=input_path,
                                force=True,
                                reference_genome=reference_genome,
                                contig_recoding=recode,
                                skip_invalid_loci=True)
                  .rows()
                  .repartition(100)
                  .key_by('locus', 'alleles')
                  )

    clinvar_tb = clinvar_tb.checkpoint(output=output_path, overwrite=overwrite)

    if export_tsv:
        (clinvar_tb
         .flatten()
         .export(f"{output_path}.tsv.bgz")
         )

    return clinvar_tb


def create_gevir_tb(
        input_path: str,
        output_path: str,
        fields: list = None,
        overwrite: bool = False,
        export_tsv: bool = False,
) -> hl.Table:
    """
    Create a Hail Table with gene-level constraint metrics from GeVir.

    :param input_path: Path to the input file.
    :param output_path: Path to save the output Hail Table.
    :param fields: List of fields to select from the input table.
    :param overwrite: Whether to overwrite the existing output file.
    :param export_tsv: Whether to export the table as a TSV file.
    :return: Hail Table
    """
    gevir_tb = hl.import_table(paths=input_path,
                               impute=True,
                               min_partitions=100,
                               key='gene_id')

    if fields is not None:
        gevir_tb = gevir_tb.select(*fields)

    gevir_tb = gevir_tb.checkpoint(output=output_path, overwrite=overwrite)

    if export_tsv:
        gevir_tb.export(output_path + '.tsv.bgz')

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
    Create a Hail Table with gene-level metrics from Ensembl.

    Table was created using the Ensembl Biomart database (GRCh38.p13) with the following query:
    <read xml file here>

    :param input_path: Path to the input file.
    :param output_path: Path to save the output Hail Table.
    :param fields: List of fields to select from the input table.
    :param canonical: Whether to filter canonical transcripts. Default: True.
    :param overwrite: Whether to overwrite the existing output file.
    :param export_tsv: Whether to export the table as a TSV file.

    :return: Hail Table
    """

    gene_tb = (hl.import_table(paths=input_path,
                               min_partitions=50)
               )

    # replace field names
    gene_tb = gene_tb.rename(ENSEMBL_BIOMART_FIELDS)

    # filter canonical transcripts
    if canonical:
        gene_tb = gene_tb.filter(gene_tb.canonical == "1", keep=True)

    # group by gene_id
    gene_tb = (gene_tb
               .group_by(gene_tb.gene_id)
               .aggregate(transcript_id=hl.agg.collect_as_set(gene_tb.transcript_id).filter(lambda x: x != ''),
                          protein_id=hl.agg.collect_as_set(gene_tb.protein_id).filter(lambda x: x != ''),
                          gene_synonym=hl.agg.collect_as_set(gene_tb.gene_synonym).filter(lambda x: x != ''),
                          gene_name=hl.agg.collect(gene_tb.gene_name).first(),
                          chromosome=hl.agg.collect(gene_tb.chromosome).first(),
                          gene_start=hl.agg.collect(gene_tb.gene_start).first(),
                          gene_end=hl.agg.collect(gene_tb.gene_end).first(),
                          gene_type=hl.agg.collect(gene_tb.gene_type).first()
                          )
               .key_by('gene_id')
               )

    if fields is not None:
        gene_tb = gene_tb.select(*fields)

    gene_tb = gene_tb.checkpoint(output=output_path, overwrite=overwrite)

    if export_tsv:
        gene_tb.export(output_path + '.tsv.bgz')

    return gene_tb

