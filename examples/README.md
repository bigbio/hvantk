# hvantk Examples

Runnable scripts and templates demonstrating hvantk workflows. For full documentation, see the [docs site](https://bigbio.github.io/hvantk/examples/).

## Prerequisites

```bash
poetry install
eval "$(poetry env activate)"
```

## Examples

| Directory | Description | Needs | Quick Start |
|-----------|-------------|-------|-------------|
| [`hgc/`](hgc/) | Joint genotyping QC and scalability benchmarks | Hail/Spark; `--extras hgc` for plots | `python examples/hgc/qc/hgc_qc_example.py` |
| [`psroc/`](psroc/) | Pathogenicity score ROC analysis with synthetic data | Hail/Spark; `--extras psroc` | `python examples/psroc/run_psroc_example.py` |
| [`enrichex/`](enrichex/) | Gene set overlap enrichment and burden testing | Overlap: base install, `--extras enrichex` for plots/reports. Burden: Hail/Spark | `python examples/enrichex/overlap_enrichment_example.py` |
| [`clinvar/`](clinvar/) | ClinVar variant streaming: filter by gene set and pathogenicity, derive training labels | Hail/Spark; network (or a local VCF with `--skip-download`) | `python examples/clinvar/clinvar_streamer_example.py` |
| [`clingen/`](clingen/) | ClinGen gene-disease validity queries and MONDO ontology categorization | Hail/Spark; network | `python examples/clingen/run_with_real_data.py` |
| [`ancestry/`](ancestry/) | Ancestry inference with PCA + Random Forest | Hail/Spark; `--extras ancestry` | `python examples/ancestry/run_ancestry_example.py` |
| [`1k_genome/`](1k_genome/) | Build 1000 Genomes NYGC reference MatrixTable | Hail/Spark; large download | `python examples/1k_genome/build_1kg_nygc.py --help` |
| [`rerank/`](rerank/) | Cohort-aware re-ranking of a burden result into a calibrated per-gene table | `--extras ml`; no Hail, no network | `hvantk rerank -c examples/rerank/config.yaml -o out/scores.tsv` |

See each subdirectory's README for prerequisites, run instructions, and expected outputs.
