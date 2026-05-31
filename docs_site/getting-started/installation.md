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
| `ml` | scikit-learn | ML-backed analyses |
| `ancestry` | scikit-learn, matplotlib, seaborn | `hvantk ancestry-inference` |
| `psroc` | scikit-learn, matplotlib, plotly | `hvantk psroc` |
| `hgc` | matplotlib, seaborn, plotly, gnomad | `hvantk hgc` QC plots/reports |
| `ptm` | cptac | CPTAC PTM downloads |
| `constraint` | tspex, matplotlib, seaborn | `hvantk ptm constraint` |
| `duckdb` | duckdb | DuckDB-backed catalog queries |

`scikit-learn` is optional — install `ml`, `ancestry`, or `psroc` if you use
those analyses.

## Prerequisites

- **Python** ≥ 3.10
- **Hail** ≥ 0.2.137
- **Java** 8 or 11 (required by Hail/Spark)
