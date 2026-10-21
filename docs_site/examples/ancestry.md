# Ancestry inference examples

Predict genetic ancestry for samples using PCA projection and Random Forest classification against a labeled reference panel (e.g., 1000 Genomes).

See the [Ancestry reference](../tools/ancestry.md) for full CLI options, the method, and caveats.

## CLI usage

```bash
# Basic inference
hvantk ancestry-inference \
  --query-mt cohort.mt \
  --reference-mt 1kg_reference.mt \
  -o results/predictions.ht \
  --output-dir results/

# With custom parameters
hvantk ancestry-inference \
  --query-mt cohort.mt \
  --reference-mt 1kg_reference.mt \
  -o results/predictions.ht \
  --output-dir results/ \
  --n-pcs 20 \
  --ancestry-col pop \
  --generate-report
```

`-o`/`--output-ht` is required; `--output-dir` only controls where the *additional*
outputs (TSV, report, model, loadings) land.

## Python API

```python
import hail as hl
from hvantk.algorithms.ancestry import run_ancestry_inference

result = run_ancestry_inference(
    query_mt=hl.read_matrix_table("cohort.mt"),
    reference_mt=hl.read_matrix_table("1kg_reference.mt"),
    ancestry_col="pop",
    n_pcs=20,
)
result.generate_report("results/ancestry_report.html")
result.save("results/ancestry")
```

`run_ancestry_inference` takes `hl.MatrixTable` objects, not paths — read them with
`hl.read_matrix_table` first. There is no `output_dir`/`generate_report` keyword on the
function itself; call `result.generate_report(...)` and `result.save(...)` separately.

## Synthetic data for testing

Generate test data using Hail's `balding_nichols_model()`. Its default column key is
`sample_idx` (an integer), so re-key to a string `s` column before use — the pipeline
assumes a column field literally named `s`:

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
combined_mt = combined_mt.annotate_cols(
    s=hl.str("sample_") + hl.str(combined_mt.sample_idx),
    ancestry=pop_labels[combined_mt.pop],
)
combined_mt = combined_mt.key_cols_by("s")

# Split into reference (labeled) and query (ancestry column dropped)
ref_mt = combined_mt.filter_cols(combined_mt.sample_idx < 200)
query_mt = combined_mt.filter_cols(combined_mt.sample_idx >= 200).drop("ancestry")
```

See [`examples/ancestry/run_ancestry_example.py`](https://github.com/bigbio/hvantk/blob/main/examples/ancestry/run_ancestry_example.py)
for the full end-to-end version of this, including validation against the labels
hidden from the query cohort.

## Expected outputs

- `predictions.ht` / `predictions.tsv` — `s`, `predicted_ancestry`, `ancestry_probability`,
  `prob_<POP>` per population, `_ancestry_source`, `_known_ancestry`, and `PC1...PCn`
  (see [Outputs](../tools/ancestry.md#outputs) for the full schema).
- `ancestry_report.html` — interactive HTML report.
- `pipeline_stats.json` — pipeline execution statistics.
- `plots/` — PCA scatter plots, ancestry proportions, confusion matrix, variance explained.

## Data requirements

- **Query MT**: Hail MatrixTable with genotype calls keyed by `(locus, alleles)` on
  rows and `s` (sample ID, `str`) on columns.
- **Reference MT**: same keying, plus a categorical ancestry-label column (name it via
  `ancestry_col`).

## Runnable scripts

See the [`examples/ancestry/`](https://github.com/bigbio/hvantk/tree/main/examples/ancestry/) directory for end-to-end examples with synthetic data generation.
