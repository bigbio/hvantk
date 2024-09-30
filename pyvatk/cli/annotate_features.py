# eam
# 06.04.22


import argparse
import sys

import hail as hl

from pyvatk.utils.annotate import (annotate_ccr,
                                  annotate_gevir,
                                  annotate_rnaseq_expression,
                                  annotate_ppi,
                                  annotate_gnomad_af,
                                  annotate_ensembl_gene,
                                  annotate_dbnsfp_scores,
                                  annotate_gnomad_constraint_metrics,
                                  annotate_degs,
                                  annotate_variant_id,
                                  annotate_clinvar_clnsig,
                                  annotate_hca)


project_dir = None
out_path = f'{project_dir}/data/features'

"""
Annotate variant table with features from multiple sources.

"""


def check_variant_tb(t: hl.Table,
                     gene_col: str):

    expected_fields = ['locus', 'alleles', gene_col]
    if all([f in t.row for f in expected_fields]):
        pass
    else:
        sys.exit(f'Expected table with at least `locus`, `alleles` and {gene_col} fields.')


def main(args):
    # Init Hail
    hl.init(default_reference='GRCh38')

    ht = hl.read_table(
        args.variant_ht
    )
    print(ht.row)

    gene_col = args.gene_col

    # check minimal requirements for input variant table.
    check_variant_tb(ht,
                     gene_col)

    # filter to bi-allelic variants
    ht = ht.filter(hl.len(ht.alleles) == 2)

    # annotate clinvar significance
    ht = annotate_clinvar_clnsig(ht)

    # annotate variant ID from locus and alleles
    ht = annotate_variant_id(ht)

    # annotate gene/transcript ensembl IDs
    ht = annotate_ensembl_gene(ht, gene_symbol_col=gene_col)

    # annotate ccr
    ht = annotate_ccr(ht)

    # annotate gvir
    ht = annotate_gevir(ht,
                        gene_id_col='GeneID')

    # annotate rnaseq expression
    ht = annotate_rnaseq_expression(ht,
                                    gene_id_col='GeneID')

    # annotate gnomad af
    ht = annotate_gnomad_af(ht)

    # annotate gnomad constraint metrics
    ht = annotate_gnomad_constraint_metrics(ht,
                                            transcript_id_col='TranscriptID')

    # annotate interactome sites
    ht = annotate_ppi(ht)

    # annotate HCA
    # ht = annotate_degs(ht,
    #                   gene_symbol_col=gene_col)
    ht = annotate_hca(ht,
                      gene_id_col='GeneID')

    # annotate deleterious scores
    ht = annotate_dbnsfp_scores(ht,
                                transcript_id_col='TranscriptID')

    # write as HT
    output_ht_path = f'{args.output_ht}/ts.denovo.features.ht'
    ht = (ht
          .checkpoint(output=output_ht_path,
                      overwrite=True)
          )

    if args.write_to_file:
        (ht
         .flatten()
         .export(f'{output_ht_path}.tsv.bgz')
         )

    # Stop Hail
    hl.stop()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--variant_ht',
                        help='Path to HailTable with variants and gene symbol keyed by `locus` and `alleles`...',
                        type=str, default=None)

    parser.add_argument('--gene_col',
                        help='Name of gene symbol column in input HailTable',
                        type=str, default=None)

    parser.add_argument('-o', '--output_ht',
                        help='Path to output HailTable with features annotations',
                        type=str, default=out_path)

    parser.add_argument('-wf', '--write_to_file', help='Write output to BGZ-compressed file',
                        action='store_true')

    args = parser.parse_args()

    main(args)

