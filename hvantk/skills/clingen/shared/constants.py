"""Plugin-local constants for the clingen skill.

Moved from hvantk/core/constants.py per issue #120 (clean-core principle).
"""

# URLs
CLINGEN_DOWNLOADS_URL = "https://search.clinicalgenome.org/kb/downloads"
CLINGEN_BASE_URL = "https://search.clinicalgenome.org/kb/gene-validity/download"
CLINGEN_FILE_PREFIX = "Clingen-Gene-Disease-Summary"
CLINGEN_HEADER_SKIP_LINES = 6

# Field-rename map: source CSV column -> snake_case
CLINGEN_GENE_DISEASE_FIELDS = {
    "GENE SYMBOL": "gene_symbol",
    "GENE ID (HGNC)": "hgnc_id",
    "DISEASE LABEL": "disease_label",
    "DISEASE ID (MONDO)": "mondo_id",
    "MOI": "mode_of_inheritance",
    "SOP": "sop_version",
    "CLASSIFICATION": "classification",
    "ONLINE REPORT": "report_url",
    "CLASSIFICATION DATE": "classification_date",
    "GCEP": "gene_curation_expert_panel",
}

# Classification levels (strongest first)
CLINGEN_CLASSIFICATION_LEVELS = [
    "Definitive",
    "Strong",
    "Moderate",
    "Limited",
    "Disputed",
    "Refuted",
    "No Known Disease Relationship",
]
