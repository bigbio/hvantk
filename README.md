[![Python Package using Conda](https://github.com/bigbio/hvantk/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/bigbio/hvantk/actions/workflows/python-package-conda.yml)
[![Python application](https://github.com/bigbio/hvantk/actions/workflows/python-app.yml/badge.svg)](https://github.com/bigbio/hvantk/actions/workflows/python-app.yml)

# hvantk

Hail-based multiomics variant annotation toolkit.

## Description

`hvankt` is a annotation toolkit that uses hail to annotate variants and genes with multiple omics data types (e.g.,
variant prediction scores, gene or protein expression). The library is designed to be modular and extensible,
allowing users to add new data types and sources. The main goal is to leverages multiomics integration and annotations
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
hvantk ucsc-matrix -e hvantk/tests/testdata/raw/ucsc/exprMatrix.test.tsv.bgz -m hvantk/tests/testdata/raw/ucsc/meta.test.tsv -o data/ucsc/exprMatrix.mt
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

### 3. Create annotation tables from raw sources:

```bash
hvantk mktables --raw_data_path /path/to/raw_data --clinvar --interactome --gevir --gnomad_metrics
```

This command creates annotation tables from raw data sources for ClinVar, interactome, GeVIR, and gnomAD metrics. Make sure to replace `/path/to/raw_data` with the actual path to your raw data directory. See [README.sources.md](README.sources.md) for instructions on how to download the raw data.

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
