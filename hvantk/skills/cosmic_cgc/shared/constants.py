"""Plugin-local constants for the cosmic_cgc skill.

Moved from hvantk/core/constants.py per issue #120 (clean-core principle).
"""

COSMIC_CGC_FILE_PREFIX = "cosmic-cgc"

# Field-rename maps, one per export header generation: source column ->
# snake_case. Each map lists its generation's whole header: the builder renames
# with both (`COSMIC_CGC_FIELDS`), then requires every column of the generation
# the header is closest to and stamps that generation's schema ID
# (`COSMIC_CGC_HEADER_GENERATIONS`; see SKILL.md s4). `build_rename_map()`
# matches case- and separator-insensitively, so "Chr Band" and "CHR_BAND" both
# match either spelling.
#
# Current (v103+) header, in export order; confirmed against a licensed
# Cosmic_CancerGeneCensus_v103_GRCh38 export (see SKILL.md s4).
COSMIC_CGC_CURRENT_FIELDS = {
    "GENE_SYMBOL": "gene_symbol",
    "NAME": "gene_name",
    "COSMIC_GENE_ID": "cosmic_gene_id",
    "CHROMOSOME": "chromosome",
    "GENOME_START": "genome_start",
    "GENOME_STOP": "genome_stop",
    "CHR_BAND": "chr_band",
    "SOMATIC": "somatic",
    "GERMLINE": "germline",
    "TUMOUR_TYPES_SOMATIC": "tumour_types_somatic",
    "TUMOUR_TYPES_GERMLINE": "tumour_types_germline",
    "CANCER_SYNDROME": "cancer_syndrome",
    "TISSUE_TYPE": "tissue_type",
    "MOLECULAR_GENETICS": "molecular_genetics",
    "ROLE_IN_CANCER": "role_in_cancer",
    "MUTATION_TYPES": "mutation_types",
    "TRANSLOCATION_PARTNER": "translocation_partner",
    "OTHER_GERMLINE_MUT": "other_germline_mut",
    "OTHER_SYNDROME": "other_syndrome",
    "TIER": "classification",
    "SYNONYMS": "synonyms",
}

# Legacy header: Title-Case names, with an Entrez ID, a combined
# "Genome Location" string and a Hallmark flag where the current header has a
# COSMIC gene ID and split coordinates.
COSMIC_CGC_LEGACY_FIELDS = {
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

# Schema ID stamped for each header generation. plugin.yaml must accept every
# key here (its schema_id plus schema_ids).
COSMIC_CGC_HEADER_GENERATIONS = {
    "cosmic-cgc-v2": COSMIC_CGC_CURRENT_FIELDS,
    "cosmic-cgc-legacy-v1": COSMIC_CGC_LEGACY_FIELDS,
}

# Every known column, for renaming before the generation is known. Built from the
# registry, so a newly registered generation's columns are renamed too.
COSMIC_CGC_FIELDS = {
    raw: field
    for fields in COSMIC_CGC_HEADER_GENERATIONS.values()
    for raw, field in fields.items()
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
