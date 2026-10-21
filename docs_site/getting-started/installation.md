# Installation

## Requirements

- **Python** 3.10, 3.11 or 3.12 (the versions CI tests).
- **Java 11.** Hail runs on the JVM and is built for Java 11; it warns on other versions.
- **[Poetry](https://python-poetry.org/) 2.x** for a development install. A plain `pip`
  install works too.

Hail (0.2.137 or later) and its bundled Spark are installed as ordinary Python
dependencies. You do not need a Spark installation: on a laptop or a single node, Hail
starts a local Spark session the first time a command needs it.

## Install with Poetry (recommended)

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
poetry install
eval "$(poetry env activate)"
```

`poetry install` uses the committed `poetry.lock`, so you get the dependency set that CI
and the container image were built from.

## Install with pip

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
pip install -e .
```

## Check the installation

```bash
hvantk utils check-install
```

It checks for a proxy setting that would break Hail, starts Hail, builds a small
synthetic MatrixTable and prints its schema. `python -m hvantk` works wherever the
`hvantk` script is not on your `PATH`.

## Optional features (extras)

The base install covers dataset builds and most analyses. Plotting, machine learning
and a few provider-specific dependencies are extras. Install the ones you need, or all
of them:

```bash
poetry install --extras "ancestry psroc"     # or: poetry install --all-extras
pip install -e ".[ancestry,psroc]"
```

| Extra | Pulls in | Needed for |
|---|---|---|
| `ancestry` | scikit-learn, matplotlib, seaborn | `hvantk ancestry-inference` |
| `psroc` | scikit-learn, matplotlib, plotly | `hvantk psroc` |
| `ml` | scikit-learn | `hvantk rerank` |
| `hgc` | matplotlib, seaborn | `hvantk hgc` QC plots and the HTML report |
| `enrichex` | matplotlib, seaborn | `hvantk enrichex` plots and `--generate-report` |
| `viz` | matplotlib, seaborn, plotly | `hvantk qtlcascade` plots and report; `hvantk ptm report` and `--save-plots` |
| `constraint` | tspex, matplotlib, seaborn, statsmodels | `hvantk ptm constraint` and `hvantk ptm test` |
| `ptm` | cptac, sorted-nearest | `hvantk download cptac-phospho` and `cptac-expression` |
| `expression` | scanpy | `hvantk expression summarize` and `markers` |
| `cohort` | | nothing: scipy is a base dependency; kept so `hvantk[cohort]` still installs |
| `interactive` | plotly | nothing: no command imports plotly today |
| `duckdb` | duckdb | nothing: no command imports duckdb today |

Without the extra, most of these commands stop with a message that names the extra to
install.

A few details for the `expression` extra:

- `hvantk expression describe` and `summarize-ucsc` work on a base install.
- `hvantk ptm constraint --expression-source anndata --expression-metric mean` also needs
  it. The default metric, `median`, does not.
- **It does not install on Intel macOS.** scanpy depends on numba and llvmlite, and
  llvmlite ships no x86_64 macOS wheel from 0.47 on. It installs normally on Linux and
  Apple Silicon. On an Intel Mac, run the expression commands on Linux or the cluster;
  the rest of the toolkit is unaffected.

## Container image for clusters

For Slurm clusters, hvantk ships an Apptainer definition
([`containers/hvantk.def`](https://github.com/bigbio/hvantk/blob/main/containers/hvantk.def))
built from `poetry.lock` with Java 11. [Running on HPC](../guide/hpc-migration.md)
covers building it and running Hail jobs with it.

## Troubleshooting

- **Hail hangs or fails with `JSONDecodeError`.** An HTTP proxy is intercepting Hail's
  localhost traffic. Unset `HTTP_PROXY`/`HTTPS_PROXY`, or set
  `NO_PROXY=localhost,127.0.0.1`. `hvantk utils check-install` warns about this.
- **Java version warnings or JVM start-up errors.** Point `JAVA_HOME` at a Java 11 JDK.
- **Out of memory in local mode.** Give the Spark driver more memory before starting
  Python: `export PYSPARK_SUBMIT_ARGS="--driver-memory 8g pyspark-shell"`.
- **Exported VCF or TSV files contain comma decimals (`AF=0,25`).** Hail formats numbers
  with the JVM's locale. Set
  `export JAVA_TOOL_OPTIONS="-Duser.language=en -Duser.country=US"`.
- **Spark log lines and warnings on stderr** (loopback hostname, native-hadoop library,
  illegal reflective access) are harmless.
