import hail as hl

from pyvatk.constants import RAW_RESOURCE_PATHS, DATA_PATH

output_dir_default = f'{DATA_PATH}/data/ht'
raw_resource_paths = RAW_RESOURCE_PATHS


def update_gnomad_constraint_metrics_tb() -> hl.Table:
    gnomad_path = raw_resource_paths.get('gnomad_metrics_path')
    gnomad_tb = hl.import_table(paths=gnomad_path,
                                impute=True,
                                min_partitions=100,
                                key='transcript')
    gnomad_tb = (gnomad_tb
                 .select(loeuf=gnomad_tb.oe_lof_upper,
                         moeuf=gnomad_tb.oe_mis_upper)
                 )
    return gnomad_tb


def update_interactome_tb() -> hl.Table:
    interactome_bed = raw_resource_paths.get('interactome_path')
    ppi_tb = (hl.import_bed(path=interactome_bed,
                            skip_invalid_intervals=True,
                            reference_genome='GRCh38',
                            )
              .repartition(100)
              .distinct()
              )
    return ppi_tb


def update_gene_ensembl_ann_tb() -> hl.Table:
    gene_ann_path = raw_resource_paths.get('gene_ann_path')
    gene_tb = (hl.import_table(paths=gene_ann_path,
                               impute=True,
                               min_partitions=100)
               )
    gene_tb = (gene_tb
               .group_by(gene_tb.GeneID,
                         gene_tb.TranscriptID,
                         gene_tb.Gene)
               .aggregate(Gene_Synonym=hl.agg.collect_as_set(gene_tb.Gene_Synonym))
               .key_by('GeneID')
               )

    return gene_tb


def update_clinvar_tb() -> hl.Table:
    clinvar_path = raw_resource_paths.get('clinvar_path')
    recode = {f"{i}": f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
    clinvar_tb = (hl.import_vcf(path=clinvar_path,
                                force=True,
                                reference_genome='GRCh38',
                                contig_recoding=recode,
                                skip_invalid_loci=True)
                  .rows()
                  .repartition(100)
                  .key_by('locus', 'alleles')
                  )

    return clinvar_tb


def update_gevir_tb() -> hl.Table:
    gevir_path = raw_resource_paths.get('gevir_path')

    gevir_tb = (hl.import_table(paths=gevir_path,
                                impute=True,
                                min_partitions=100,
                                key='gene_id')
                )

    return gevir_tb


def update_scell_deg_tb() -> hl.Table:
    scell_path = raw_resource_paths.get('scell_heart_path')

    sdeg_tb = (hl.import_table(paths=scell_path,
                               impute=True,
                               min_partitions=50)
               )
    sdeg_tb = (sdeg_tb
               .group_by('gene')
               .aggregate(cluster_id=hl.agg.collect_as_set(sdeg_tb.cluster_id),
                          cluster_name=hl.agg.collect_as_set(sdeg_tb.Cluster_name))
               .key_by('gene')
               )
    return sdeg_tb


def update_hca_tb() -> hl.Table:
    hca_path = raw_resource_paths.get('scell_hca_path')

    hca_tb = hl.import_table(paths=hca_path,
                             min_partitions=100,
                             impute=True)

    cell_categories_rank = ["adipocyte",
                            "atrial_cardiomyocyte",
                            "endothelial",
                            "fibroblast",
                            "lymphoid",
                            "mesothelial",
                            "myeloid",
                            "neuronal",
                            "not_assigned",
                            "pericyte",
                            "smooth_muscle_cell",
                            "ventricular_cardiomyocyte",
                            "doublet"]

    hca_tb = (hca_tb
              .transmute(gene_name=hca_tb.name,
                         gene_id=hca_tb.name2.split("[.]")[0],
                         gene_id_version=hca_tb.name2.split("[.]")[1],
                         expScores=hca_tb.expScores.split("[,]").map(lambda x: hl.float(x))
                         )
              .select('gene_id', 'expScores')
              )

    hca_tb = hca_tb.transmute(mean_umi=hl.dict(hl.zip(cell_categories_rank, hca_tb.expScores)))

    ks = hca_tb.mean_umi.key_set().collect()[0]
    hca_tb = (hca_tb
              .transmute(**{k: hca_tb.mean_umi.get(k)
                            for k in ks})
              .key_by('gene_id')
              )

    hca_tb.show()
    hca_tb.describe()

    return hca_tb


def update_rnaseq_tb() -> hl.Table:
    rna_table = raw_resource_paths.get('rnaseq_path')

    tb = hl.import_table(paths=rna_table,
                         min_partitions=100,
                         delimiter='\\s',
                         quote='"')

    # getting all expression level (per sample) fields
    expr_fields = [f for f in tb.row if f != 'Gene']

    # parse expression value
    tb = tb.annotate(
        **{f: hl.parse_float(tb[f]) for f in expr_fields}
    ).key_by('Gene')

    # convert Table -> MatrixTable
    mt = tb.to_matrix_table_row_major(columns=expr_fields,
                                      entry_field_name='cpm',
                                      col_field_name='sample')

    # add sample meta-data
    mt_ann = mt.annotate_cols(_sample_array=mt.sample.split("[.]"))
    mt_ann = mt_ann.annotate_cols(organ=mt_ann._sample_array[0],
                                  time_point=mt_ann._sample_array[1],
                                  sample_index=mt_ann._sample_array[2],
                                  organismus='Human')

    # annotate developmental stages info grouping by time points
    dev = [f'{t}wpc' for t in range(1, 9)]
    mat = [f'{t}wpc' for t in range(9, 25)]

    dev_stage_ann_expr = {'dev_stage':
                              hl.case()
                              .when(hl.set(dev).contains(mt_ann.time_point), 'development')
                              .when(hl.set(mat).contains(mt_ann.time_point), 'maturation')
                              .default('postnatal')
                          }
    mt_ann = mt_ann.annotate_cols(**dev_stage_ann_expr)

    # annotate mean expression per time points and developmental stages
    mt_t = (mt_ann
            .annotate_rows(mean_expr_time_point=hl.agg.group_by(hl.struct(organ=mt_ann.organ,
                                                                          time_point=mt_ann.time_point),
                                                                hl.agg.mean(mt_ann.cpm)),
                           mean_expr_dev_stage=hl.agg.group_by(hl.struct(organ=mt_ann.organ,
                                                                         dev_stage=mt_ann.dev_stage),
                                                               hl.agg.mean(mt_ann.cpm)))

            )

    return mt_t.rows()
