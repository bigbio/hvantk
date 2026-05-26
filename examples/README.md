# hvantk Examples

Runnable scripts and templates demonstrating hvantk workflows. For full documentation, see the [docs site](https://bigbio.github.io/hvantk/examples/).

## Prerequisites

```bash
poetry install
eval "$(poetry env activate)"
```

## Examples

| Directory | Description | Quick Start |
|-----------|-------------|-------------|
| [`hgc/`](hgc/) | Joint genotyping QC and scalability benchmarks | `python examples/hgc/qc/hgc_qc_example.py` |
| [`psroc/`](psroc/) | Pathogenicity score ROC analysis with synthetic data | `python examples/psroc/run_psroc_example.py` |
| [`enrichex/`](enrichex/) | Gene set overlap enrichment and burden testing | `python examples/enrichex/overlap_enrichment_example.py` |
| [`clinvar/`](clinvar/) | ClinVar data streaming and chunk processing | `python examples/clinvar/clinvar_streamer_example.py` |
| [`clingen/`](clingen/) | ClinGen gene-disease validity queries and ontology categorization | `python examples/clingen/run_with_real_data.py` |
| [`ancestry/`](ancestry/) | Ancestry inference with PCA + Random Forest | `python examples/ancestry/run_ancestry_example.py` |
| [`1k_genome/`](1k_genome/) | Build 1000 Genomes NYGC reference MatrixTable | `python examples/1k_genome/build_1kg_nygc.py --help` |

See each subdirectory's README for prerequisites, run instructions, and expected outputs.
