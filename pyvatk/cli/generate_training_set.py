# eam
# 07.04.22

"""

Generate training set from Clinvar

"""

import hail as hl

from fschd.utils.data_utils import (get_chd_gene_set,
                                    get_clinvar_ht)
from settings import project_dir

out_dir = f"{project_dir}/data/training_set"

print(out_dir)

hl.init()

# Pathogenic labels from Clinvar
PATHOGENIC_LABEL_CLINVAR = ['Pathogenic/Likely_pathogenic',
                            'Likely_pathogenic',
                            'Pathogenic']

# Disease-specific labels from Clinvar
CHD_LABEL_CLINVAR = ['Congenital_heart_disease',
                     'Congenital_heart_defect']

# Benign labels from Clinvar
BENIGN_LABEL_CLINVAR = ['Benign/Likely_benign',
                        'Likely_benign',
                        'Benign']

# import known CHD-associated genes
chd_gene_set = get_chd_gene_set()

# import clinvar data base
clinvar_ht = get_clinvar_ht()

# add gene column to clinvar table
clinvar_ht = (clinvar_ht
              .annotate(gene=clinvar_ht.info.GENEINFO.split("[:]")[0])
              )

# filter to non-synonymous variants
clinvar_ht = (clinvar_ht
              .annotate(Consequence=
                        clinvar_ht.info.MC.map(lambda x:
                                               x.split("[|]")[1])
                        )
              )

clinvar_ht = (clinvar_ht
              .filter(~clinvar_ht.Consequence.any(lambda x: x == "synonymous_variant"))
              )

# annotate TP/FP label based on Pathogenic/Benign annotations from clinvar
ts_ann_expr = {'is_tp_site':
                   hl.case()
                       .when(clinvar_ht.info.CLNSIG.any(lambda x: hl.set(PATHOGENIC_LABEL_CLINVAR).contains(x)) &
                             hl.set(chd_gene_set).contains(clinvar_ht.gene), True)
                       .when(clinvar_ht.info.CLNDN.any(lambda x: hl.set(CHD_LABEL_CLINVAR).contains(x)), True)
                       .default(False),
               'is_tn_site':
                   hl.case()
                       .when(clinvar_ht.info.CLNSIG.any(lambda x: hl.set(BENIGN_LABEL_CLINVAR).contains(x)), True)
                       .default(False)
               }

ts_ht = (clinvar_ht
         .annotate(**ts_ann_expr)
         )

ts_ht = (ts_ht
         .filter(ts_ht.is_tp_site != ts_ht.is_tn_site)
         )

ts_ht = (ts_ht
         .annotate(rf_label=hl.case()
                   .when(ts_ht.is_tp_site, 'TP')
                   .when(ts_ht.is_tn_site, 'TN')
                   .or_missing())
         )

ts_ht = ts_ht.filter(hl.is_defined(ts_ht.rf_label))

ts_ht = (ts_ht
         .select('gene', 'rf_label')
         )

# export results
ht_out_path = f"{out_dir}/ts.clinvar.ht"
ts_ht = (ts_ht
         .checkpoint(output=ht_out_path,
                     overwrite=True)
         )

ts_ht.export(
    f"{ht_out_path}.tsv"
)

hl.stop()
