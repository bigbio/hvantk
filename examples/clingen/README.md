# ClinGen Examples

Examples for querying ClinGen Gene-Disease Validity data and categorizing diseases by ontology.

For full documentation, see the [ClinGen examples guide](https://bigbio.github.io/hvantk/examples/clingen/).

## Contents

- **`run_with_real_data.py`** - Read the ClinGen Hail Table built by the Prerequisites command below (exits 1 if it is missing), extract genesets per disease, view classification stats
- **`run_ontology_categorization.py`** - Categorize diseases using MONDO ontology hierarchy

## Prerequisites

Run these from the **repository root**:

```bash
# Download today's ClinGen export and build the Hail Table both scripts read
hvantk reprocess clingen:gene-disease \
  --raw-dir examples/clingen/data \
  --output examples/clingen/results/clingen_full.ht
```

`run_ontology_categorization.py` downloads the MONDO ontology (`data/mondo.obo`) on its
first run.

## Quick Start

```bash
python examples/clingen/run_with_real_data.py
python examples/clingen/run_ontology_categorization.py
```

## Expected Outputs

**`run_with_real_data.py`** (in `results/`):
- `clingen_full.ht/` - Hail Table with ClinGen data
- `genesets_per_disease_full.json` - Genes grouped by disease label
- `category_summary_full.tsv` - Summary of genes per disease category
- `classification_distribution_full.tsv` - Distribution across classification levels

**`run_ontology_categorization.py`** (in `results/`):
- `ontology_category_summary.tsv` - Summary of genes per MONDO category
- `genes_by_ontology_category.json` - Genes grouped by ontology category
- `diseases_by_ontology_category.json` - Diseases grouped by ontology category

## Requirements

- Hail installed and configured
- Downloaded ClinGen CSV and MONDO OBO files (see Prerequisites)
