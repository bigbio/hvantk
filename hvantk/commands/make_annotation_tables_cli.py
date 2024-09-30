"""
Generate annotation tables from multiple (raw) sources

"""

import click

from hvantk.settings import CONTEXT_SETTINGS, RAW_DATA_PATH, set_raw_data_path
from hvantk.utils.make_tables import (create_interactome_tb,
                                      create_rnaseq_tb,
                                      create_clinvar_tb,
                                      create_gevir_tb,
                                      create_scell_deg_tb,
                                      create_hca_tb,
                                      create_gene_ensembl_ann_tb,
                                      create_gnomad_constraint_gene_metrics_tb)

output_dir_default = f'{RAW_DATA_PATH}/annotation_tables'


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """A package for gene and variant annotation."""
    pass


def make_annotation_tables_from_raw_sources(raw_data_path: str,
                                            ccr: bool = False,  # TODO: implement CCR table
                                            interactome: bool = False,
                                            temporal_rnaseq: bool = False,
                                            clinvar: bool = False,
                                            gevir: bool = False,
                                            scell_heart_deg: bool = False,
                                            hca_rnaseq: bool = False,
                                            gene_ensembl: bool = False,
                                            gnomad_metrics: bool = False,
                                            output_dir: str = output_dir_default,
                                            default_ref_genome: str = 'GRCh38'):
    # set the raw data path
    set_raw_data_path(raw_data_path)

    if interactome:
        bed_ppi = create_interactome_tb()
        bed_ppi.checkpoint(
            f'{output_dir}/interactome.{default_ref_genome}.ht',
            overwrite=True
        )

    if temporal_rnaseq:
        rnaseq_tb = create_rnaseq_tb()
        rnaseq_tb.checkpoint(
            f'{output_dir}/rnaseq.human.ht',
            overwrite=True
        )

    if clinvar:
        clinvar_tb = create_clinvar_tb()
        clinvar_tb.checkpoint(
            f'{output_dir}/clinvar.{default_ref_genome}.ht',
            overwrite=True
        )

    if gevir:
        gevir_tb = create_gevir_tb()
        gevir_tb.checkpoint(
            f'{output_dir}/gevir.metrics.ht',
            overwrite=True
        )

    if scell_heart_deg:
        deg_tb = create_scell_deg_tb()
        deg_tb.checkpoint(
            f'{output_dir}/scell.heart.degs.ht',
            overwrite=True
        )

    if hca_rnaseq:
        hca_tb = create_hca_tb()
        hca_tb.checkpoint(
            f'{output_dir}/hca.heart.ht',
            overwrite=True
        )

    if gene_ensembl:
        gene_tb = create_gene_ensembl_ann_tb()
        gene_tb.checkpoint(
            f'{output_dir}/gene.ann.ensembl.ht',
            overwrite=True
        )

    if gnomad_metrics:
        gnomad_tb = create_gnomad_constraint_gene_metrics_tb()
        gnomad_tb.checkpoint(
            f'{output_dir}/gnomad.metrics.ht',
            overwrite=True
        )


@click.command('mktables', short_help='Create annotation tables from raw sources.')
@click.option('--raw_data_path', default=RAW_DATA_PATH, type=str, required=True,
              help='Path to raw data directory.')
@click.option('--output_dir', default=output_dir_default, type=str, required=True,
              help='Output directory to copy created Hail tables')
@click.option('--ccr',
              is_flag=True, help='Create/update CCR table from source.')
@click.option('--interactome',
              is_flag=True, help='Create/update CCR table from source.')
@click.option('--temporal_rnaseq',
              is_flag=True, help='Create/update RNAseq table from source.')
@click.option('--clinvar',
              is_flag=True, help='Create/update Clinvar table from source.')
@click.option('--gevir',
              is_flag=True, help='Create/update GeVIR score table from raw source.')
@click.option('--scell_heart_deg',
              is_flag=True, help='Create/update table with DEGs from cardiac-specific cell clusters')
@click.option('--hca_rnaseq',
              is_flag=True, help='Create/update table with gene/cell expression levels from HCA dataset (UCSC)')
@click.option('--gene_ensembl',
              is_flag=True, help='Create/update gene annotation table from Ensembl.')
@click.option('--gnomad_metrics',
              is_flag=True, help='Create/update transcript-specific constraint metrics from gnomad database')
@click.option('--default_ref_genome', default='GRCh38', type=str,
              help='Default reference genome to start Hail. Only GRCh38 is supported for now.')
@click.pass_context
def make_annotation_tables_cli(ctx, raw_data_path,  output_dir, ccr, interactome, temporal_rnaseq, clinvar, gevir,
                               scell_heart_deg, hca_rnaseq, gene_ensembl, gnomad_metrics, default_ref_genome):

    # exit if no flat parameter is set
    if not any([ccr, interactome, temporal_rnaseq, clinvar,
                gevir, scell_heart_deg, hca_rnaseq, gene_ensembl, gnomad_metrics]):
        click.echo('No flag set. Please set at least one flag to create/update a table.')
        ctx.abort()

    make_annotation_tables_from_raw_sources(raw_data_path,
                                            ccr,
                                            interactome,
                                            temporal_rnaseq,
                                            clinvar,
                                            gevir,
                                            scell_heart_deg,
                                            hca_rnaseq,
                                            gene_ensembl,
                                            gnomad_metrics,
                                            output_dir,
                                            default_ref_genome)


if __name__ == '__main__':
    cli()
