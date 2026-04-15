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

# ---------------------------------------------------------------------------
# Phase-2 defaults (notebook A / M / K / N)
# ---------------------------------------------------------------------------

# Proximal flanking window for SYMBOL-based variant annotation:
# ±7 aa in codon space ≈ ±21 bp in genomic coordinates (notebook N Cell 11).
# One-sided; applied to both codon_start and codon_end (see ptm.annotate).
PROXIMAL_BP = 21

# Default MAF thresholds referenced by Phase-2 documentation.
# Notebooks M and K do NOT apply MAF cut-offs (they work from ClinVar B/LB
# variants only). Notebook N uses MAF <= 0.001 on an internal CHD cohort.
# These values are kept as Phase-2 defaults/placeholders for any
# caller that wants a consistent set of bins (common/rare/ultra-rare);
# callers that need the exact notebook-N rare filter should pass
# ``rare=1e-3`` explicitly.
DEFAULT_MAF_THRESHOLDS = {
    "common": 0.01,
    "rare": 0.001,
    "ultra_rare": 1e-5,
}

# Ordered category for the binned-interaction LMM (notebook K).
# "b0_none" is the reference (zero-expression) bin; Q1-Q4 are quartiles of
# log2(expr + 1) across positive-expression variants. Actual bin count in
# a fit may be smaller when pd.qcut collapses ties (duplicates="drop").
EXPRESSION_BIN_LABELS = ["b0_none", "b1_Q1", "b2_Q2", "b3_Q3", "b4_Q4"]

# Pseudocount for log10(AF + eps) transformation used by notebooks K and M.
LOG_AF_EPSILON = 1e-8

# Per-stratum filter thresholds for the constraint LMM (notebook M Cell 5).
LMM_MIN_N_PTM = 30
LMM_MIN_N_NONPTM = 30
LMM_MIN_MIXED_GENES = 10

# Sparsity gates for the binned-interaction LMM (notebook K Cell 4d).
LMM_BINNED_MIN_POS_EXPR = 100
LMM_BINNED_MIN_CELL_N = 5
