"""Constants and default values for PTM module."""

# UniProt REST API
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_API_FIELDS = "accession,gene_names,ft_mod_res,xref_ensembl,sequence"
UNIPROT_HUMAN_PTM_QUERY = "organism_id:9606 AND reviewed:true AND ft_mod_res:*"
UNIPROT_BATCH_SIZE = 500

# Ensembl GTF
ENSEMBL_RELEASE = "113"
ENSEMBL_GTF_URL = (
    f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}"
    f"/gtf/homo_sapiens/Homo_sapiens.GRCh38.{ENSEMBL_RELEASE}.gtf.gz"
)
ENSEMBL_GTF_FILENAME = f"Homo_sapiens.GRCh38.{ENSEMBL_RELEASE}.gtf.gz"

# PeptideAtlas Phospho Build
PEPTIDEATLAS_PHOSPHO_BASE_URL = "https://peptideatlas.org/builds/human/phospho"
PEPTIDEATLAS_LATEST_BUILD_DATE = "202512"
PEPTIDEATLAS_LATEST_BUILD_ID = "606"

# CPTAC Phospho (via cptac Python package)
CPTAC_CANCER_TYPES = [
    "brca", "ccrcc", "coad", "gbm",
    "hnscc", "lscc", "luad", "ov", "pdac", "ucec",
]
CPTAC_CANCER_CLASS_MAP = {
    "brca": "Brca",
    "ccrcc": "Ccrcc",
    "coad": "Coad",
    "gbm": "Gbm",
    "hnscc": "Hnscc",
    "lscc": "Lscc",
    "luad": "Luad",
    "ov": "Ov",
    "pdac": "Pdac",
    "ucec": "Ucec",
}

# PTM type categories (UniProt MOD_RES description prefixes)
PTM_TYPE_CATEGORIES = {
    "Phosphoserine": "phosphorylation",
    "Phosphothreonine": "phosphorylation",
    "Phosphotyrosine": "phosphorylation",
    "N6-acetyllysine": "acetylation",
    "N6-methyllysine": "methylation",
    "N6,N6-dimethyllysine": "methylation",
    "N6,N6,N6-trimethyllysine": "methylation",
    "Asymmetric dimethylarginine": "methylation",
    "Symmetric dimethylarginine": "methylation",
    "Omega-N-methylarginine": "methylation",
    "N-linked (GlcNAc...) asparagine": "n_glycosylation",
    "O-linked (GalNAc...) threonine": "o_glycosylation",
    "O-linked (GalNAc...) serine": "o_glycosylation",
    "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO)": "sumoylation",
    "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO2)": "sumoylation",
    "Glycyl lysine isopeptide (Lys-Gly)": "ubiquitination",
}

# Default flanking window sizes to test (in residues)
DEFAULT_FLANKING_WINDOWS = [3, 5, 7, 10]
DEFAULT_FLANKING_CODONS = 5

# Output TSV column names
PTM_OUTPUT_COLUMNS = [
    "chrom",
    "codon_start",
    "codon_end",
    "strand",
    "uniprot_id",
    "gene_symbol",
    "residue_pos",
    "amino_acid",
    "ptm_type",
    "ptm_category",
    "source_db",
    "evidence_type",
    "n_observations",
]

# Transcript resolution method names (for logging/QC)
RESOLUTION_METHODS = ["xref_mane", "xref_any", "gene_mane"]
