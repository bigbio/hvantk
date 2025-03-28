# Constants for the project

# Ensembl biomart fields mapping
# ['Gene stable ID', 'Transcript stable ID', 'Protein stable ID', 'Chromosome/scaffold name',
# 'Gene start (bp)', 'Gene end (bp)', 'Ensembl Canonical', 'Gene name', 'Gene type', 'Gene Synonym']

ENSEMBL_BIOMART_FIELDS = {
    "Gene stable ID": "gene_id",
    "Transcript stable ID": "transcript_id",
    "Protein stable ID": "protein_id",
    "Chromosome/scaffold name": "chromosome",
    "Gene start (bp)": "gene_start",
    "Gene end (bp)": "gene_end",
    "Ensembl Canonical": "canonical",
    "Gene name": "gene_name",
    "Gene type": "gene_type",
    "Gene Synonym": "gene_synonym",
}

# UCSC Cell Browser base URL
UCSC_CELL_BROWSER_BASE_URL = "https://cells.ucsc.edu"
EXPRESSION_MATRIX_FILE_NAME = "exprMatrix.tsv.gz"
METADATA_FILE_NAME = "meta.tsv"
