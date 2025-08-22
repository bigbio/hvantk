[![Python Package using Conda](https://github.com/bigbio/hvantk/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/bigbio/hvantk/actions/workflows/python-package-conda.yml)
[![Python application](https://github.com/bigbio/hvantk/actions/workflows/python-app.yml/badge.svg)](https://github.com/bigbio/hvantk/actions/workflows/python-app.yml)

# hvantk

Hail-based multiomics variant annotation toolkit.

## Description

`hvantk` is an annotation toolkit that uses Hail to annotate variants and genes with multiple omics data types (e.g.,
variant prediction scores, gene or protein expression). The library is designed to be modular and extensible,
allowing users to add new data types and sources. The main goal is to leverage multiomics integration and annotations
from heterogeneous sources to improve the interpretation of genetic variants.

## Installation

Download the source code and install the package using Poetry:

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
poetry install
```

If you don't have Poetry installed, you can install it using:

```bash
pip install poetry
```

or, if you prefer conda:

```bash
conda install -c conda-forge poetry
```

Then, activate the environment:

```bash
poetry shell
```

## Usage Examples

### 1. Download UCSC Cell Browser data:

```bash
hvantk ucsc-downloader --dataset adultPancreas --output-dir data/ucsc
```

This command downloads the expression matrix and metadata for the `adultPancreas` dataset from the UCSC Cell Browser and saves it to the `data/ucsc` directory.

### 2. Convert UCSC Cell data to Hail matrix table:

```bash
hvantk mkmatrix ucsc -e hvantk/tests/testdata/raw/ucsc/exprMatrix.test.tsv.bgz -m hvantk/tests/testdata/raw/ucsc/meta.test.tsv -o data/ucsc/exprMatrix.mt
```

This command converts the expression matrix and metadata files from the UCSC Cell Browser into a Hail matrix table format.

Example matrix table schema:

```markdown
---

Global fields:
None

---

Column fields:
'cell_id': str
'metadata': struct {
orig_ident: str,
nCount_RNA: int32,
nFeature_RNA: int32,
percent_mt: float64,
Rep: int32,
Age: int32,
Region: str,
RNA_snn_res_0_8: int32,
seurat_clusters: int32,
clusters: int32,
colors: str,
major_cell_class: str
}

---

Row fields:
'gene': str

---

Entry fields:
'x': int32

---

Column key: ['cell_id']
Row key: ['gene']

---
```

### 3A. Batch-create annotation tables from a recipe:

Create a recipe JSON (YAML also supported if PyYAML is installed):

```json
{
  "tables": [
    {
      "name": "clinvar",
      "input": "/data/clinvar_2024.vcf.bgz",
      "output": "/out/clinvar.ht",
      "params": {"reference_genome": "GRCh38", "export_tsv": true}
    },
    {
      "name": "interactome",
      "input": "/data/insider.bed.bgz",
      "output": "/out/interactome.ht",
      "params": {"reference_genome": "GRCh38"}
    }
  ]
}
```

Run:

```bash
hvantk mktable-batch --recipe /path/to/tables.json
```

### 3B. Create a single table from an explicit raw file:

- ClinVar (VCF → HT keyed by locus, alleles):

```bash
hvantk mktable clinvar --raw-input /path/to/clinvar.vcf.bgz --output-ht /path/to/clinvar.ht --ref-genome GRCh38 --overwrite
```

- Interactome (BED intervals → HT keyed by interval):

```bash
hvantk mktable interactome --raw-input /path/to/interactome.bed.bgz --output-ht /path/to/interactome.ht
```

- GeVIR (TSV keyed by gene_id):

```bash
hvantk mktable gevir --raw-input /path/to/gevir.tsv.bgz --output-ht /path/to/gevir.ht --fields oe_syn_upper,oe_mis_upper
```

- gnomAD constraint metrics (TSV keyed by gene_id):

```bash
hvantk mktable gnomad-metrics --raw-input /path/to/gnomad.tsv.bgz --output-ht /path/to/gnomad.ht
```

- Ensembl gene annotations (Biomart TSV keyed by gene_id):

```bash
hvantk mktable ensembl-gene --raw-input /path/to/biomart.tsv.bgz --output-ht /path/to/ensembl.ht --no-canonical
```

Run `hvantk mktable --help` or `hvantk mktable <subcommand> --help` for full options.

### 2B. Batch-create MatrixTables from a recipe:

Create a recipe JSON (YAML also supported if PyYAML is installed):

```json
{
  "matrices": [
    {
      "name": "ucsc",
      "inputs": {
        "expression_matrix": "/data/ucsc/expr.tsv.bgz",
        "metadata": "/data/ucsc/meta.tsv"
      },
      "output": "/out/ucsc.mt",
      "params": {"gene_column": "gene", "overwrite": true}
    }
  ]
}
```

Run:

```bash
hvantk mkmatrix-batch --recipe /path/to/matrices.json
```

For more examples and recipes, see docs/USAGE.md and examples/recipes/.

## Annotation sources

A full description of the sources and how to download the data is available in the
[README.sources.md](README.sources.md) file.

- Variants and genomic regions

  - Missense variants prediction scores (from dbNSFP)
  - ClinVar annotations
  - gnomAD annotations (e.g. allele frequencies)
  - Protein-protein interaction site (INSIDER)
  - Ensemble gene annotations
  - GeVIR score (PMID:31873297)
  - Coding-constrained region (CCR) score

- Bulk RNA-seq data

  - Human tissue expression (brain, heart, liver, kidney), multiple developmental time points (E-MTAB-6814)

- Single-cell RNA-seq data

  - Embryonic human heart single-cell RNA-seq data (PMID:31835037).
  - Human heart single-cell RNA-seq data (PMID:31835037).
  - Human heart cell atlas (UCSC, https://doi.org/10.1038/s41586-020-2797-4).

- Protein expression data
  - TODO: Add protein expression data sources.

# Things to do:

- Add a section to download the data from the sources.
- Add a section about conversion from local files. including local mapping files of they are needed.
- Some small benchmarks with loom -> to the annotation tool in hail.

## Developer quickstart

- Install and activate the environment (see Installation), then run tests:
  - `pytest -q`
- Explore the CLI to see available commands:
  - `hvantk --help`
- Typical workflow when adding a new data source:
  1) Define a data product contract (Table/MatrixTable schema + metadata)
  2) Write a downloader (optional) and a builder that outputs a Hail Table/MatrixTable
  3) Register the dataset in a small manifest (provenance, versions, hashes)
  4) Create streamers (transformers) and compose a recipe to answer a biological question
  5) Add tiny tests using the fixtures in hvantk/tests/testdata

See:
- docs/DEVELOPING.md – dev workflow and contracts
- docs/STREAMERS_AND_RECIPES.md – streamer interface and JSON recipes (YAML optional)
- docs/DATA_CATALOG.md – hosting strategy and dataset registry format

## Limitations and strategy

- Heterogeneous omics, full-table builds: Prefer a slice-first approach. Builders should support selectors (genes/regions, tissues/cell types, timepoints) so users don’t have to build everything. Cache by parameter hash to reuse slices.
- Limited hosting: Use a lightweight data catalog (JSON first; YAML optional) that points to immutable remote URIs (S3/GCS/Zenodo/DOI) with checksums. Host only manifests and small indices in this repo.
- Streamers and pipelines: Define a tiny plugin contract for streamers (read -> transform -> write) and compose them with JSON "recipes" (YAML optional). Keep streamers stateless and testable on tiny fixtures.

Quick starts for each topic and examples live under examples/ (see examples/recipes/ and examples/datasets/).
