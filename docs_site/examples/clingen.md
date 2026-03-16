# ClinGen Examples

Examples demonstrating how to work with ClinGen Gene-Disease Validity data using hvantk.

## Overview

ClinGen (Clinical Genome Resource) provides curated gene-disease validity assessments. The hvantk toolkit provides tools to:

- Build Hail Tables from ClinGen CSV data
- Query genes by disease, classification, or mode of inheritance
- Categorize diseases using MONDO ontology hierarchy
- Generate summary statistics and reports
- Export gene sets for downstream analysis

## Quick Start

```bash
# Download ClinGen data and MONDO ontology
cd examples/clingen
mkdir -p data

# Download ClinGen Gene-Disease Validity CSV
curl -L --fail "https://search.clinicalgenome.org/kb/gene-validity/download" -o data/clingen_gene_disease.csv

# Download MONDO ontology (for ontology-based categorization)
curl -L --fail "https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.obo" -o data/mondo.obo

# Run example with real data
python run_with_real_data.py

# Run example with ontology-based disease categorization
python run_ontology_categorization.py
```

## Examples

### `run_with_real_data.py`

Basic example demonstrating:
1. Building a ClinGen Hail Table from downloaded CSV data
2. Extracting genesets grouped by individual diseases
3. Generating summaries of genes per disease category (keyword-based)
4. Viewing classification distribution and statistics

**Output files (in `results/`):**
- `clingen_full.ht/` - Hail Table with ClinGen data
- `genesets_per_disease_full.json` - Genes grouped by disease label
- `category_summary_full.tsv` - Summary of genes per disease category
- `classification_distribution_full.tsv` - Distribution across classification levels

### `run_ontology_categorization.py`

Advanced example demonstrating ontology-based disease categorization using MONDO:
1. Building a ClinGen Hail Table from downloaded CSV data
2. Loading and parsing the MONDO disease ontology
3. Categorizing diseases based on their ontological hierarchy (is_a relationships)
4. Comparing keyword-based vs ontology-based categorization approaches

**Output files (in `results/`):**
- `ontology_category_summary.tsv` - Summary of genes per MONDO category
- `genes_by_ontology_category.json` - Genes grouped by ontology category
- `diseases_by_ontology_category.json` - Diseases grouped by ontology category

## ClinGen Classification Levels

From highest to lowest confidence:
1. **Definitive** - Conclusive evidence for gene-disease relationship
2. **Strong** - Strong evidence supporting the relationship
3. **Moderate** - Moderate evidence supporting the relationship
4. **Limited** - Limited evidence, requires more research
5. **Disputed** - Evidence is conflicting
6. **Refuted** - Evidence refutes the relationship
7. **No Known Disease Relationship**

## Using the ClinGenStreamer API

```python
from hvantk.data.clingen_streamer import ClinGenStreamer

# Initialize streamer with a built Hail Table
streamer = ClinGenStreamer("path/to/clingen.ht", init_hail=False)
streamer.setup()

# Get genes by classification level
definitive_genes = streamer.get_genes_by_classification(
    min_classification="Definitive",
    as_set=True
)

# Get genes by disease term (keyword matching)
cancer_genes = streamer.get_genes_by_disease(
    disease_terms=["cancer", "carcinoma"],
    match_mode="contains",
    min_classification="Moderate",
    as_set=True
)

# Get genesets per disease
genesets = streamer.get_geneset_per_disease(min_classification="Moderate")

# Ontology-based categorization (requires MONDO OBO file)
results = streamer.categorize_by_ontology(
    mondo_obo_path="data/mondo.obo",
    min_classification="Limited"
)
# results["cardiovascular disease"]["genes"] -> set of genes

# Get summary DataFrame of ontology categories
summary_df = streamer.categorize_by_ontology_summary(
    mondo_obo_path="data/mondo.obo"
)

# Get summary statistics
stats = streamer.compute_stats()
print(f"Unique genes: {stats['unique_genes']}")
print(f"Unique diseases: {stats['unique_diseases']}")
```

## MONDO Ontology Categories

The ontology-based categorization uses the MONDO disease ontology hierarchy.
Default categories include:

- **Organ/System-based**: cardiovascular, nervous system, metabolic, immune system, eye, respiratory, etc.
- **Cancer**: neoplasm, cancer
- **Genetic**: hereditary disease, autosomal dominant/recessive, X-linked
- **Developmental**: developmental disorder, neurodevelopmental disorder, malformation syndrome

See `hvantk/utils/mondo_parser.py` for the full list of `MONDO_DISEASE_CATEGORIES`.

## Data Sources

- **ClinGen**: https://search.clinicalgenome.org/kb/gene-validity
- **MONDO Ontology**: https://github.com/monarch-initiative/mondo

## Requirements

- Hail installed and configured
- hvantk installed (`poetry install`)
- Downloaded ClinGen CSV and MONDO OBO files (see Quick Start)
