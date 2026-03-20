# Data Sources

This page covers all annotation and expression data sources supported by hvantk: what they are, where to get them, and how to build Hail Tables from the raw data. Sources are split into two categories: those with **built-in downloaders** (automated) and those that require **manual download** (too large, license-gated, or fragile URLs).

For building Hail Tables and MatrixTables from downloaded data, see the [Usage Guide](usage.md).

## File format note

Downloaded `.gz` files may be standard gzip (single-threaded in Hail) rather than BGZF (parallel). Use `--auto-convert-bgz` during table or matrix builds, or pre-convert with:

```bash
hvantk convert-bgz input.gz
```

## Sources with built-in downloaders

| Source | Command | Approx. Size |
|---|---|---|
| ClinVar | `hvantk clinvar-downloader` | ~500 MB |
| ClinGen | `hvantk clingen-downloader` | ~5 MB |
| HGNC | `hvantk hgnc-downloader` | ~20 MB |
| UCSC Cell Browser | `hvantk ucsc-downloader` | varies |
| Expression Atlas | `hvantk expression-atlas-downloader` | varies |

### ClinVar

Clinically relevant variants and their annotations (e.g. Pathogenic, Benign, VUS).
URL: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/

```bash
# Download latest ClinVar VCF (GRCh38) with tabix index
hvantk clinvar-downloader --output-dir data/clinvar

# Download a specific archived version
hvantk clinvar-downloader --version 20260101 --output-dir data/clinvar

# Download GRCh37 build, verify checksum
hvantk clinvar-downloader --genome-build GRCh37 --verify-md5
```

### ClinGen

```bash
# Download today's ClinGen Gene-Disease Validity snapshot
# Output: Clingen-Gene-Disease-Summary-<YYYY-MM-DD>.csv
hvantk clingen-downloader --output-dir data/clingen

# Check download availability
hvantk clingen-downloader --list-versions
```

### HGNC

```bash
# Download HGNC complete gene nomenclature set
hvantk hgnc-downloader --output-dir data/hgnc
```

### UCSC Cell Browser

The UCSC Cell Browser hosts 267+ datasets. About half are **collections** (groups of
related datasets with no expression matrix at the top level). Use `--list_datasets`
and `--search` to discover downloadable datasets.

```bash
# Discover available datasets
hvantk ucsc-downloader --list_datasets

# Search by name, organism, or tissue (expands collections to show children)
hvantk ucsc-downloader --list_datasets --search heart
hvantk ucsc-downloader --list_datasets --search pancreas

# Download a leaf dataset directly
hvantk ucsc-downloader --dataset adultPancreas --output-dir data/ucsc

# Download a child dataset from a collection (use the full path)
hvantk ucsc-downloader --dataset hoc/all-heart --output-dir data/ucsc
```

> **Note:** Collection names (e.g., `hoc`) cannot be downloaded directly — they
> contain no expression matrix. Use `--search` to find child dataset paths like
> `hoc/all-heart`, then download those.

### Expression Atlas

```bash
# Download bulk RNA-seq experiments
hvantk expression-atlas-downloader --download_path data/expression_atlas
```

## Manual download sources

These sources are too large, require license acceptance, or have complex download procedures. Follow the instructions below, then use `hvantk mktable` to build Hail Tables.

### dbNSFP (~45 GB)

A database of functional prediction scores for human missense variants.
URL: https://sites.google.com/site/jpopgen/dbNSFP

**Download**: Requires academic license acceptance. Download from the project page:
https://sites.google.com/site/jpopgen/dbNSFP

**Pre-processing**: dbNSFP is distributed as per-chromosome `.gz` files (standard gzip, not BGZF). The builder expects a **single combined file**, so concatenate and BGZF-compress first:

```bash
# Concatenate per-chromosome files into a single BGZF file
# (header is taken from chr1; remaining files skip the header line)
head -1 <(zcat dbNSFP4.9a_variant.chr1.gz) > /tmp/dbnsfp_header.txt
(cat /tmp/dbnsfp_header.txt && for f in dbNSFP4.9a_variant.chr*.gz; do zcat "$f" | tail -n +2; done) \
  | bgzip -@ 4 > dbNSFP4.9a_variant.bgz
```

**Build**:

```bash
# Option 1: Pre-converted BGZF (recommended)
hvantk mktable dbnsfp \
  --raw-input dbNSFP4.9a_variant.bgz \
  --output-ht dbnsfp.ht

# Option 2: Auto-convert during build (requires a single already-merged .gz)
# This only works if you have already concatenated the per-chromosome files
# into a single gzip file (e.g., dbNSFP4.9a_variant.gz).
# --auto-convert-bgz re-compresses the single .gz as BGZF; it does NOT
# assemble per-chromosome files.
hvantk mktable dbnsfp \
  --raw-input dbNSFP4.9a_variant.gz \
  --output-ht dbnsfp.ht \
  --auto-convert-bgz
```

### gnomAD constraint metrics (~50 MB for gene-level)

Gene-level constraint metrics (pLI, LOEUF, missense Z-score) from gnomAD v4.1.
URL: https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv

**Download**:

```bash
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
```

**Build**:

```bash
hvantk mktable gnomad-metrics \
  --raw-input gnomad.v4.1.constraint_metrics.tsv \
  --output-ht gnomad_metrics.ht
```

### INSIDER interactome (~100 MB)

Protein-protein interaction sites from the INSIDER database.
URL: http://interactomeinsider.yulab.org/downloads.html

**Download**: Visit http://interactomeinsider.yulab.org/downloads.html and download the interaction site BED file.

**Build**:

```bash
hvantk mktable interactome \
  --raw-input insider_interaction_sites.bed.bgz \
  --output-ht interactome.ht
```

### Ensembl gene annotations (~800 MB)

Ensembl gene annotations (gene name, gene ID, biotype, transcript ID).
URL: https://www.ensembl.org/info/data/ftp/index.html

**Download**: Export from BioMart with the required attributes matching `ENSEMBL_BIOMART_FIELDS` in `hvantk/core/constants.py`. Alternatively, download from the Ensembl FTP:
https://www.ensembl.org/info/data/ftp/index.html

**Build**:

```bash
hvantk mktable ensembl-gene \
  --raw-input biomart_export.tsv.bgz \
  --output-ht ensembl_gene.ht
```

### GeVIR (~20 GB)

Gene variation intolerance ranking scores.
URL: https://www.nature.com/articles/s41588-019-0560-2

**Download**: Supplementary data from the Nature publication:
https://www.nature.com/articles/s41588-019-0560-2

**Build**:

```bash
hvantk mktable gevir \
  --raw-input gevir_metrics.tsv.bgz \
  --output-ht gevir.ht
```

### CCR - Coding-Constrained Regions (~50 MB)

Highly constrained coding regions in the human genome.
URL: https://www.nature.com/articles/s41588-018-0294-6

**Download**: Supplementary data from the Nature publication above.

**Note**: No builder is currently available for CCR. This is planned for a future release.

## Expression data sources

These datasets are used to build expression MatrixTables via the UCSC Cell Browser and Expression Atlas downloaders.

### Bulk RNA-seq

- **Human tissue expression E-MTAB-6814** — Human tissue gene expression (brain, heart, liver, kidney), multiple developmental time points.
  URL: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6814

### Single-cell RNA-seq

- **Human heart scRNA-seq (Asp 2019)** — Embryonic human heart single-cell RNA-seq data 6.5 wpc (PMID:31835037).
  URL: https://data.mendeley.com/datasets/mbvhhf8m62/2
- **Human heart scRNA-seq (Farah 2024)** — Single-cell RNA-seq data of the developing human heart, 9-15 wpc.
  URL: https://cells.ucsc.edu/?bp=heart&ds=hoc
- **Human heart cell atlas (HCA)** — Adult human heart cell atlas (https://doi.org/10.1038/s41586-020-2797-4).
  URL: https://cells.ucsc.edu/?bp=heart&ds=heart-cell-atlas
