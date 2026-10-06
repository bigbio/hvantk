"""Plugin-local constants for the cosmic_cgc skill.

Moved from hvantk/core/constants.py per issue #120 (clean-core principle).
"""

COSMIC_CGC_FILE_PREFIX = "cosmic-cgc"

# Field-rename map: source CSV column -> snake_case
#
# Covers two header generations. `build_rename_map()` matches case- and
# separator-insensitively, so e.g. both "Chr Band" (legacy) and "CHR_BAND"
# (current) normalize to the same key and need only one entry below.
#
# Current (v103+) exports additionally carry COSMIC_GENE_ID, CHROMOSOME,
# GENOME_START, GENOME_STOP, which have no legacy-header equivalent (the
# legacy export instead carried a single combined "Genome Location" string
# and a separate "Entrez GeneId" -- neither of which the current export
# includes). Both generations' entries are kept so either header maps
# correctly; confirmed against a licensed Cosmic_CancerGeneCensus_v103_GRCh38
# export (see SKILL.md s4).
COSMIC_CGC_FIELDS = {
    "Gene Symbol": "gene_symbol",
    "Name": "gene_name",
    "COSMIC_GENE_ID": "cosmic_gene_id",
    "CHROMOSOME": "chromosome",
    "GENOME_START": "genome_start",
    "GENOME_STOP": "genome_stop",
    # Legacy-export-only fields below this line (absent from v103+; see
    # SKILL.md s4/s5) -- EXCEPT "Tier" onward, which are shared by both
    # header generations (normalized matching handles the case/punctuation
    # difference, e.g. "Chr Band" vs "CHR_BAND").
    "Entrez GeneId": "entrez_id",
    "Genome Location": "genome_location",
    "Hallmark": "hallmark",
    # Shared by both header generations from here on.
    "Tier": "classification",
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
