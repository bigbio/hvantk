# Core constants for the project

from pathlib import Path
import logging

logger = logging.getLogger(__name__)

# Base directory of the package (hvantk)
BASE_DIR = Path(__file__).resolve().parent
logger.debug(f"Base directory: {BASE_DIR}")

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
    "HGNC_DOWNLOAD_URL",
    "HGNC_INFO_URL",
    "HGNC_GENE_FIELDS",
    "HGNC_PIPE_SEPARATED_FIELDS",
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
