# eam
# 29.03.22

# path (global) variable settings

import os


# context settings for click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

# global project path
# DATA_PATH = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = None

# raw resources paths
RAW_RESOURCE_PATHS = {
   'interactome_path':      f'{DATA_PATH}/resources/interactome/Interactome_INSIDER_hg38_stripped.bed',
   'clinvar_path':          f'{DATA_PATH}/resources/clinvar/clinvar_20220403.vcf.gz',
   'rnaseq_path':           f'{DATA_PATH}/resources/rnaseq-expression/E-MTAB-6814.Human.CPM.txt',
   'gene_ann_path':         f'{DATA_PATH}/resources/ensembl/gene.ensembl.canonical.042022.tsv',
   'gnomad_metrics_path':   f'{DATA_PATH}/resources/gnomad/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz',
   'gevir_path':            f'{DATA_PATH}/resources/gevir/gevir_metrics_pmid31873297.tsv.txt',
   'scell_heart_path':      f'{DATA_PATH}/resources/rnaseq-expression/deg_scell_heart_pmid31835037.tsv',
   'scell_hca_path':        f'{DATA_PATH}/resources/rnaseq-expression/hca_cells_ucsc_042022.tsv'
}


# processed resources paths
PROCESSED_RESOURCE_PATHS = {}

