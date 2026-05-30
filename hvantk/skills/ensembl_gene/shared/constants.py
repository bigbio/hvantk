"""Plugin-local constants for the ensembl_gene skill.

Moved from hvantk/core/constants.py per issue #120 (clean-core principle).
"""

# Ensembl/biomart field-rename map: BioMart export column name -> snake_case
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
