"""Plugin-local constants for the cosmic_cgc skill.

Moved from hvantk/core/constants.py per issue #120 (clean-core principle).
"""

COSMIC_CGC_FILE_PREFIX = "cosmic-cgc"

# Field-rename map: source CSV column -> snake_case
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
