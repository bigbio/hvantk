# Core constants for the project

from pathlib import Path
import logging

logger = logging.getLogger(__name__)

# Base directory of the package (hvantk)
BASE_DIR = Path(__file__).resolve().parent
logger.debug(f"Base directory: {BASE_DIR}")

## UCSC Cell Browser base URL
UCSC_CELL_BROWSER_BASE_URL = "https://cells.ucsc.edu"
logger.debug(f"UCSC Cell Browser base URL: {UCSC_CELL_BROWSER_BASE_URL}")
EXPRESSION_MATRIX_FILE_NAME = "exprMatrix.tsv.gz"
logger.debug(f"Expression matrix file name: {EXPRESSION_MATRIX_FILE_NAME}")
METADATA_FILE_NAME = "meta.tsv"
logger.debug(f"Metadata file name: {METADATA_FILE_NAME}")

# Path to the JSON file containing the UCSC cell datasets
UCSC_JSON_FILE_PATH = BASE_DIR.parent / "resources" / "cells_ucsc_datasets.json"
logger.debug(f"UCSC JSON file path: {UCSC_JSON_FILE_PATH}")

# UCSC gene and cell ID columns
UCSC_CELL_ID_COLUMN = "cell_id"
UCSC_GENE_COLUMN = "gene"

## Expression Atlas base URL
EXPRESSION_ATLAS_BASE_URL = "https://www.ebi.ac.uk/gxa/experiments-content"
logger.debug(f"Expression Atlas base URL: {EXPRESSION_ATLAS_BASE_URL}")

# Path to the new unified registry system
REGISTRY_ROOT_PATH = BASE_DIR.parent / "resources" / "registry"

# Backward compatibility - path to transcriptomics datasets (replaces expression_atlas.json)
TRANSCRIPTOMICS_DATASETS_PATH = REGISTRY_ROOT_PATH / "transcriptomics" / "datasets.json"

# Legacy path for compatibility (deprecated)
EXPRESSION_ATLAS_JSON_FILE_PATH = (
    BASE_DIR.parent / "resources" / "expression_atlas.json"
)

logger.debug(f"Registry root path: {REGISTRY_ROOT_PATH}")
logger.debug(f"Transcriptomics datasets path: {TRANSCRIPTOMICS_DATASETS_PATH}")
logger.debug(f"Legacy Expression Atlas path: {EXPRESSION_ATLAS_JSON_FILE_PATH}")

# ClinGen Gene-Disease Validity
CLINGEN_DOWNLOADS_URL = "https://search.clinicalgenome.org/kb/downloads"
CLINGEN_BASE_URL = "https://search.clinicalgenome.org/kb/gene-validity/download"
CLINGEN_FILE_PREFIX = "Clingen-Gene-Disease-Summary"
CLINGEN_HEADER_SKIP_LINES = 6

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

CLINGEN_CLASSIFICATION_LEVELS = [
    "Definitive",
    "Strong",
    "Moderate",
    "Limited",
    "Disputed",
    "Refuted",
    "No Known Disease Relationship",
]

logger.debug(f"ClinGen base URL: {CLINGEN_BASE_URL}")

# GenCC (Gene Curation Coalition)
GENCC_BASE_URL = (
    "https://thegencc.org/download/action/submissions-export-tsv?format=new"
)
GENCC_FILE_PREFIX = "gencc-submissions"

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

logger.debug(f"GenCC base URL: {GENCC_BASE_URL}")

# COSMIC Cancer Gene Census (CGC)
COSMIC_CGC_FILE_PREFIX = "cosmic-cgc"

COSMIC_CGC_FIELDS = {
    "Gene Symbol": "gene_symbol",
    "Name": "gene_name",
    "Entrez GeneId": "entrez_id",
    "Genome Location": "genome_location",
    "Tier": "classification",
    "Hallmark": "hallmark",
    "Chr Band": "chr_band",
    "Somatic": "somatic",
    "Germline": "germline",
    "Tumour Types(Somatic)": "tumour_types_somatic",
    "Tumour Types(Germline)": "tumour_types_germline",
    "Cancer Syndrome": "cancer_syndrome",
    "Tissue Type": "tissue_type",
    "Molecular Genetics": "molecular_genetics",
    "Role in Cancer": "role_in_cancer",
    "Mutation Types": "mutation_types",
    "Translocation Partner": "translocation_partner",
    "Other Germline Mut": "other_germline_mut",
    "Other Syndrome": "other_syndrome",
    "Synonyms": "synonyms",
}

COSMIC_CGC_CLASSIFICATION_LEVELS = [
    "Tier 1",
    "Tier 2",
]

COSMIC_TISSUE_TYPES = {
    "E": "Epithelial",
    "L": "Lymphoid",
    "M": "Mesenchymal",
    "O": "Other",
}

COSMIC_MUTATION_CONTEXTS = ["somatic", "germline", "both"]

# HGNC Gene Nomenclature
HGNC_DOWNLOAD_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
HGNC_INFO_URL = "https://www.genenames.org/download/statistics-and-files/"

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

logger.debug(f"HGNC download URL: {HGNC_DOWNLOAD_URL}")

# ClinVar VCF downloads
CLINVAR_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"
CLINVAR_FTP_BASE_GRCh37 = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37"

# ClinVar clinical significance labels
CLINVAR_PATHOGENIC_LABELS = [
    "Pathogenic/Likely_pathogenic",
    "Likely_pathogenic",
    "Pathogenic",
]
CLINVAR_BENIGN_LABELS = [
    "Benign/Likely_benign",
    "Likely_benign",
    "Benign",
]

logger.debug(f"ClinVar FTP base URL: {CLINVAR_FTP_BASE}")

# AlphaGenome defaults
ALPHAGENOME_DEFAULT_INTERVAL_SIZE = 1_048_576  # 1Mbp
ALPHAGENOME_DEFAULT_DENSITY_WINDOW = 50_000    # 50kb
ALPHAGENOME_DEFAULT_RETRY_BACKOFF = 2.0
ALPHAGENOME_DEFAULT_MAX_RETRIES = 3
ALPHAGENOME_DEFAULT_REQUEST_TIMEOUT = 120

# Explicit public API for this module
__all__ = [
    "BASE_DIR",
    "UCSC_CELL_BROWSER_BASE_URL",
    "EXPRESSION_MATRIX_FILE_NAME",
    "METADATA_FILE_NAME",
    "UCSC_JSON_FILE_PATH",
    "UCSC_CELL_ID_COLUMN",
    "UCSC_GENE_COLUMN",
    "EXPRESSION_ATLAS_BASE_URL",
    "REGISTRY_ROOT_PATH",
    "TRANSCRIPTOMICS_DATASETS_PATH",
    "EXPRESSION_ATLAS_JSON_FILE_PATH",  # Backward compatibility
    "CLINGEN_DOWNLOADS_URL",
    "CLINGEN_BASE_URL",
    "CLINGEN_FILE_PREFIX",
    "CLINGEN_HEADER_SKIP_LINES",
    "CLINGEN_GENE_DISEASE_FIELDS",
    "CLINGEN_CLASSIFICATION_LEVELS",
    "HGNC_DOWNLOAD_URL",
    "HGNC_INFO_URL",
    "HGNC_GENE_FIELDS",
    "HGNC_PIPE_SEPARATED_FIELDS",
    "GENCC_BASE_URL",
    "GENCC_FILE_PREFIX",
    "GENCC_SUBMISSION_FIELDS",
    "GENCC_CLASSIFICATION_LEVELS",
    "COSMIC_CGC_FILE_PREFIX",
    "COSMIC_CGC_FIELDS",
    "COSMIC_CGC_CLASSIFICATION_LEVELS",
    "COSMIC_TISSUE_TYPES",
    "COSMIC_MUTATION_CONTEXTS",
    "CLINVAR_FTP_BASE",
    "CLINVAR_FTP_BASE_GRCh37",
    "CLINVAR_PATHOGENIC_LABELS",
    "CLINVAR_BENIGN_LABELS",
    "ALPHAGENOME_DEFAULT_INTERVAL_SIZE",
    "ALPHAGENOME_DEFAULT_DENSITY_WINDOW",
    "ALPHAGENOME_DEFAULT_RETRY_BACKOFF",
    "ALPHAGENOME_DEFAULT_MAX_RETRIES",
    "ALPHAGENOME_DEFAULT_REQUEST_TIMEOUT",
]
