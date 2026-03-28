# Ancestry Inference Examples

Predict genetic ancestry for samples using PCA projection and Random Forest classification against a labeled reference panel (e.g., 1000 Genomes).

See the [Ancestry reference](../tools/ancestry.md) for full CLI options and pipeline details.

## CLI Usage

```bash
# Basic inference
hvantk ancestry-inference \
  --query-mt cohort.mt \
  --reference-mt 1kg_reference.mt \
  --output-dir results/

# With custom parameters
hvantk ancestry-inference \
  --query-mt cohort.mt \
  --reference-mt 1kg_reference.mt \
  --output-dir results/ \
  --n-pcs 20 \
  --ancestry-col pop \
  --generate-report
```

## Python API

```python
from hvantk.ancestry.pipeline import run_ancestry_pipeline

run_ancestry_pipeline(
    query_mt_path="cohort.mt",
    reference_mt_path="1kg_reference.mt",
    output_dir="results/",
    n_pcs=20,
    generate_report=True,
)
```

## Synthetic Data for Testing

Generate test data using Hail's `balding_nichols_model()`:

```python
import hail as hl

combined_mt = hl.balding_nichols_model(
    n_populations=5,
    n_samples=300,
    n_variants=10000,
    pop_dist=[0.2, 0.2, 0.2, 0.2, 0.2],
    fst=[0.12, 0.15, 0.10, 0.08, 0.11],
)

pop_labels = hl.literal(["EUR", "AFR", "EAS", "SAS", "AMR"])
combined_mt = combined_mt.annotate_cols(ancestry=pop_labels[combined_mt.pop])

# Split into reference and query
ref_mt = combined_mt.filter_cols(combined_mt.sample_idx < 200)
query_mt = combined_mt.filter_cols(combined_mt.sample_idx >= 200)
```

## Expected Outputs

- `ancestry_report.html` - Interactive HTML report
- `predictions.tsv` - Ancestry predictions (sample_id, predicted, probability)
- `pipeline_stats.json` - Pipeline execution statistics
- `plots/` - PCA scatter plots, ancestry proportions, confusion matrix, variance explained

## Data Requirements

- **Query MT**: Hail MatrixTable with genotype calls keyed by `(locus, alleles)`
- **Reference MT**: Hail MatrixTable with ancestry labels in a column annotation

## Runnable Scripts

See the [`examples/ancestry/`](https://github.com/bigbio/hvantk/tree/main/examples/ancestry/) directory for end-to-end examples with synthetic data generation.
