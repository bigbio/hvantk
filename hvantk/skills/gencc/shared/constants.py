"""Plugin-local constants for the gencc skill.

Moved from hvantk/core/constants.py per issue #120 (clean-core principle).
"""

# URLs
GENCC_BASE_URL = (
    "https://thegencc.org/download/action/submissions-export-tsv?format=new"
)
GENCC_FILE_PREFIX = "gencc-submissions"

# Field-rename map: source TSV column -> snake_case
GENCC_SUBMISSION_FIELDS = {
    "sgc_id": "sgc_id",
    "gene_curie": "hgnc_id",
    "gene_symbol": "gene_symbol",
    "disease_curie": "mondo_id",
    "disease_title": "disease_label",
    "disease_original_curie": "disease_original_id",
    "disease_original_title": "disease_original_label",
    "classification_title": "classification",
    "moi_title": "mode_of_inheritance",
    "submitter_title": "submitter",
    "submitted_as_date": "submission_date",
    "submitted_as_public_report_url": "report_url",
    "submitted_as_pmids": "pmids",
}

# Classification levels (strongest first)
GENCC_CLASSIFICATION_LEVELS = [
    "Definitive",
    "Strong",
    "Moderate",
    "Supportive",
    "Limited",
    "Disputed Evidence",
    "Refuted Evidence",
    "No Known Disease Relationship",
]
