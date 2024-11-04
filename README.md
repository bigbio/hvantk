# hvantk
Hail-based multiomics variant annotation toolkit.


## Description

`hvankt` is a annotation toolkit that uses hail to annotate variants and genes with multiple omics data types (e.g.,
variant prediction scores, gene or protein expression). The library is designed to be modular and extensible,
allowing users to add new data types and sources. The main goal is to leverages multiomics integration and annotations 
from heterogeneous sources to improve the interpretation of genetic variants.

## Installation

Download the source code and install the package using one of the following methods.

Option 1: Install via _pip_.

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
pip install . -r requirements.txt
```

Option 2: Install via _conda_.

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
conda env create -f environment.yml
conda activate hvantk-env
pip install . -r requirements.txt
```

## Annotation sources

A full description of the sources and how to download the data is available in the 
[README.sources.md](README.sources.md) file.

* Variants and genomic regions
  - Missense variants prediction scores (from dbNSFP)
  - ClinVar annotations
  - gnomAD annotations (e.g. allele frequencies)
  - Protein-protein interaction site (INSIDER)
  - Ensemble gene annotations
  - GeVIR score (PMID:31873297)
  - Coding-constrained region (CCR) score


* Bulk RNA-seq data
  - Human tissue expression (brain, heart, liver, kidney), multiple developmental time points (E-MTAB-6814)


* Single-cell RNA-seq data
  - Embryonic human heart single-cell RNA-seq data (PMID:31835037).
  - Human heart single-cell RNA-seq data (PMID:31835037).
  - Human heart cell atlas (UCSC, https://doi.org/10.1038/s41586-020-2797-4).


* Protein expression data
  - TODO: Add protein expression data sources.

# Things to do: 

- Add a section to download the data from the sources. 
- Add a section about conversion from local files. including local mapping files of they are needed. 
- Some small benchmarks with loom -> to the annotation tool in hail.

