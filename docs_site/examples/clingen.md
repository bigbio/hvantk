# ClinGen Examples

Query ClinGen Gene-Disease Validity data, extract gene sets, and categorize diseases using MONDO ontology.

## Quick Start

The ClinGen plugin has a built-in downloader, so one command downloads today's snapshot into `data/` and builds the Hail Table:

```bash
hvantk reprocess clingen:gene-disease \
  --raw-dir data/ \
  --output clingen.ht
```

If you already downloaded `Clingen-Gene-Disease-Summary-<YYYY-MM-DD>.csv` into `data/`, add `--skip-download` to reuse it.

## ClinGenStreamer API

`ClinGenStreamer` inherits from `GeneDiseaseValidityStreamer`, which provides generic query, aggregation, and integration methods. The same API is available via `GenCCStreamer` for GenCC data.

```python
from hvantk.data.clingen_streamer import ClinGenStreamer

streamer = ClinGenStreamer("clingen.ht", init_hail=False)
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

# Get genesets grouped by disease
genesets = streamer.get_geneset_per_disease(min_classification="Moderate")
```

## Ontology-Based Disease Categorization

Categorize diseases using the MONDO ontology hierarchy (is_a relationships) instead of keyword matching:

```python
# Requires MONDO OBO file
results = streamer.categorize_by_ontology(
    ontology="data/mondo.obo",
    min_classification="Limited"
)
# results["cardiovascular disease"]["genes"] -> set of genes

# Summary DataFrame
summary_df = streamer.categorize_by_ontology_summary(
    ontology="data/mondo.obo"
)
```

Default MONDO categories include organ/system-based (cardiovascular, nervous system, metabolic), cancer, genetic (hereditary, autosomal dominant/recessive), and developmental. See `hvantk/utils/mondo_parser.py` for the full list.

## Classification Levels

From highest to lowest confidence:

1. **Definitive** - Conclusive evidence
2. **Strong** - Strong evidence
3. **Moderate** - Moderate evidence
4. **Limited** - Limited evidence, requires more research
5. **Disputed** - Conflicting evidence
6. **Refuted** - Evidence refutes the relationship

## Data Sources

- **ClinGen**: <https://search.clinicalgenome.org/kb/gene-validity>
- **MONDO Ontology**: <https://github.com/monarch-initiative/mondo>

## Runnable Scripts

See the [`examples/clingen/`](https://github.com/bigbio/hvantk/tree/main/examples/clingen/) directory for complete examples with real data and ontology categorization.
