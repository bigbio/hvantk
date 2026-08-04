# Installation

## Using Poetry (recommended)

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
poetry install
eval "$(poetry env activate)"
```

## Using pip

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
pip install -e .
```

## Optional features (extras)

The core install is lightweight; visualization, machine-learning, and a few
provider-specific dependencies ship as optional extras. Install the ones you
need:

```bash
# pip
pip install -e ".[viz,hgc]"

# poetry
poetry install --extras "viz hgc"      # or: poetry install --all-extras
```

| Extra | Pulls in | Needed for |
|---|---|---|
| `viz` | matplotlib, seaborn, plotly | plots and HTML reports |
| `interactive` | plotly | interactive QC dashboards |
| `ml` | scikit-learn, scipy | ML-backed analyses |
| `ancestry` | scikit-learn, matplotlib, seaborn, scipy | `hvantk ancestry-inference` |
| `psroc` | scikit-learn, matplotlib, plotly, scipy | `hvantk psroc` |
| `hgc` | matplotlib, seaborn | `hvantk hgc` QC plots/reports |
| `ptm` | cptac, sorted-nearest | CPTAC PTM downloads |
| `constraint` | tspex, matplotlib, seaborn, scipy | `hvantk ptm constraint` |
| `enrichex` | scipy, matplotlib, seaborn | `hvantk enrichex overlap` / `burden` |
| `cohort` | scipy | `hvantk cohort burden` |
| `duckdb` | duckdb | DuckDB-backed catalog queries |
| `expression` | scanpy, scipy | `hvantk expression summarize` and `markers` |

Three command paths need this extra:

- `hvantk expression summarize` — via `summarize_expression_ad`
- `hvantk expression markers` — scanpy's `rank_genes_groups`
- `hvantk ptm constraint --expression-source anndata --expression-metric mean` —
  routes through `summarize_expression_ad` for the mean aggregation only. The
  default metric is `median`, which uses a direct numpy path and does not need
  scanpy.

`hvantk expression describe` and `summarize-ucsc` do not touch scanpy and work on
a base install. Without the extra, these paths exit with an actionable message
naming it, not a traceback.

> **`expression` is unavailable on Intel macOS.** scanpy pulls `numba` →
> `llvmlite`, which ships no x86_64 macOS wheel from 0.47 onward and fails to
> build from source. The extra installs normally on Linux and Apple Silicon.
> This is why scanpy is an extra rather than a base dependency: as a base dep it
> made the whole package uninstallable on Intel Macs. On Intel, run those two
> commands on the cluster; the rest of the toolkit is unaffected.

`scikit-learn` is optional — install `ml`, `ancestry`, or `psroc` if you use
those analyses.

## Prerequisites

- **Python** ≥ 3.10
- **Hail** ≥ 0.2.137
- **Java** 8 or 11 (required by Hail/Spark)
