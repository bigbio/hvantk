"""Plugin-local constants for the hgnc skill.

Moved from hvantk/core/constants.py per issue #120 (clean-core principle).
"""

# Download URLs
HGNC_DOWNLOAD_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
HGNC_INFO_URL = "https://www.genenames.org/download/statistics-and-files/"

# Field-rename map: source TSV column -> snake_case
HGNC_GENE_FIELDS = {
    # Core identifiers
    "hgnc_id": "hgnc_id",
    "symbol": "gene_symbol",
    "name": "gene_name",
    "status": "status",
    # Symbol history
    "alias_symbol": "alias_symbols",
    "alias_name": "alias_names",
    "prev_symbol": "prev_symbols",
    "prev_name": "prev_names",
    # Cross-references
    "ensembl_gene_id": "ensembl_gene_id",
    "entrez_id": "entrez_id",
    "uniprot_ids": "uniprot_ids",
    "refseq_accession": "refseq_id",
    "ucsc_id": "ucsc_id",
    "ccds_id": "ccds_id",
    # Classification
    "locus_group": "locus_group",
    "locus_type": "locus_type",
    "gene_group": "gene_group",
    "gene_group_id": "gene_group_id",
    # Location
    "location": "location",
    "location_sortable": "location_sortable",
    # Disease/clinical
    "omim_id": "omim_id",
    "orphanet": "orphanet_id",
    "gencc": "gencc",
    "mane_select": "mane_select",
    # Metadata
    "date_approved_reserved": "date_approved",
    "date_modified": "date_modified",
    "date_symbol_changed": "date_symbol_changed",
}

# Fields that arrive as pipe-separated multi-values and should split to array<str>
HGNC_PIPE_SEPARATED_FIELDS = [
    "alias_symbols",
    "alias_names",
    "prev_symbols",
    "prev_names",
    "uniprot_ids",
    "ccds_id",
    "gene_group",
    "gene_group_id",
    "omim_id",
    "mane_select",
]
