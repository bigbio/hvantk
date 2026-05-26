# Data Sources

> **Heads up — examples below need refresh.** This page references the retired `hvantk mktable` and `hvantk mkmatrix` CLIs. The unified replacement is `hvantk reprocess <plugin>:<dataset>` — see the [Usage Guide](usage.md#1-build-a-dataset-with-hvantk-reprocess) for the current pattern. Download commands on this page are still accurate; the build commands need adapting.

This page covers all annotation and expression data sources supported by hvantk: what they are, where to get them, and how to build datasets from the raw data. Sources are split into two categories: those with **built-in downloaders** (automated) and those that require **manual download** (too large, license-gated, or fragile URLs).

For building datasets from downloaded data, see the [Usage Guide](usage.md).

## File format note

Downloaded `.gz` files may be standard gzip (single-threaded in Hail) rather than BGZF (parallel). Use `--auto-convert-bgz` during table or matrix builds, or pre-convert with:

```bash
hvantk utils convert-bgz input.gz
```

## Sources with built-in downloaders

| Source | Command | Approx. Size |
|---|---|---|
| ClinVar | `hvantk download clinvar` | ~500 MB |
| ClinGen | `hvantk download clingen` | ~5 MB |
| GenCC | `hvantk download gencc` | ~10 MB |
| HGNC | `hvantk download hgnc` | ~20 MB |
| UCSC Cell Browser | `hvantk download ucsc` | varies |
| Expression Atlas | `hvantk download expression-atlas` | varies |

### ClinVar

Clinically relevant variants and their annotations (e.g. Pathogenic, Benign, VUS).
URL: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/

```bash
# Download latest ClinVar VCF (GRCh38) with tabix index
hvantk download clinvar --output-dir data/clinvar

# Download a specific archived version
hvantk download clinvar --version 20260101 --output-dir data/clinvar

# Download GRCh37 build, verify checksum
hvantk download clinvar --genome-build GRCh37 --verify-md5
```

### ClinGen

```bash
# Download today's ClinGen Gene-Disease Validity snapshot
# Output: Clingen-Gene-Disease-Summary-<YYYY-MM-DD>.csv
hvantk download clingen --output-dir data/clingen

# Check download availability
hvantk download clingen --list-versions
```

### GenCC

GenCC (Gene Curation Coalition) aggregates gene-disease validity assertions from 12+ submitting organizations (ClinGen, PanelApp, G2P, Orphanet, etc.).

```bash
# Download today's GenCC submissions snapshot
hvantk download gencc --output-dir data/gencc

# Check download availability
hvantk download gencc --list-versions
```

### HGNC

```bash
# Download HGNC complete gene nomenclature set
hvantk download hgnc --output-dir data/hgnc
```

**Build Hail Table**:

```bash
hvantk mktable hgnc \
  --raw-input data/hgnc/hgnc_complete_set.tsv \
  --output-ht hgnc.ht
```

### UCSC Cell Browser

The UCSC Cell Browser hosts 267+ datasets. About half are **collections** (groups of
related datasets with no expression matrix at the top level). Use `--list_datasets`
and `--search` to discover downloadable datasets.

```bash
# Discover available datasets
hvantk download ucsc --list_datasets

# Search by name, organism, or tissue (expands collections to show children)
hvantk download ucsc --list_datasets --search heart
hvantk download ucsc --list_datasets --search pancreas

# Download a leaf dataset directly
hvantk download ucsc --dataset adultPancreas --output-dir data/ucsc

# Download a child dataset from a collection (use the full path)
hvantk download ucsc --dataset hoc/all-heart --output-dir data/ucsc
```

> **Note:** Collection names (e.g., `hoc`) cannot be downloaded directly — they
> contain no expression matrix. Use `--search` to find child dataset paths like
> `hoc/all-heart`, then download those.

### Expression Atlas

```bash
# Download bulk RNA-seq experiments
hvantk download expression-atlas --download_path data/expression_atlas
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

### COSMIC Cancer Gene Census

Gene-level cancer annotations from the COSMIC Cancer Gene Census.
URL: https://cancer.sanger.ac.uk/census

**Download**: Requires COSMIC account. Download the Cancer Gene Census TSV from the COSMIC website.

**Build**:

```bash
hvantk mktable cosmic-cgc \
  --raw-input cancer_gene_census.tsv \
  --output-ht cosmic_cgc.ht
```

### UniProt PTM Sites

Curated post-translational modification sites (phosphorylation, ubiquitination, acetylation, etc.) for reviewed human proteins from UniProt/Swiss-Prot.
URL: https://www.uniprot.org/

The `hvantk ptm build` command downloads PTM data automatically via the UniProt REST API. For pre-download or manual acquisition:

**Download** (optional, for offline use):

```bash
# UniProt PTM TSV via REST API
hvantk download uniprot-ptm --output-dir data/ptm/

# Ensembl GTF for coordinate mapping (download manually)
# wget https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.*.gtf.gz -P data/ref/
```

**Build**:

```bash
# Automatic download and build
hvantk ptm build --output-dir data/ptm/ --output-ht data/ptm/ptm_sites.ht

# With pre-downloaded files
hvantk ptm build \
  --gtf-path data/ref/Homo_sapiens.GRCh38.113.gtf.gz \
  --ptm-tsv data/ptm/uniprot-ptm-human.tsv \
  --output-dir data/ptm/ \
  --output-ht data/ptm/ptm_sites.ht
```

### Ensembl GTF (~50 MB compressed)

Gene annotation with exon coordinates and CDS phases, used by the PTM mapper for residue-to-genomic coordinate mapping.
URL: https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/

**Download**: Auto-downloaded by `hvantk ptm build`. For manual download:

```bash
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
```

## QTL data

These datasets are used to build eQTL and pQTL Hail Tables for the QTL cascade pipeline.

### GTEx eQTL data

Expression quantitative trait loci from the GTEx project. Available as significant pairs (genome-wide significant associations) and allpairs (full summary statistics for coloc).

**GTEx v11** (recommended):
URL: https://www.gtexportal.org/home/downloads/adult-gtex/qtl

```bash
# Download significant pairs (Parquet format, ~50 MB per tissue)
# Navigate to GTEx Portal → Downloads → Adult GTEx → QTL → eQTL → Significant pairs

# Build significant-pairs table
hvantk mktable eqtl \
  --raw-input /data/gtex_v11/Liver.v11.signif_pairs.parquet \
  --output-ht eqtl_liver.ht \
  --source gtex_v11 \
  --tissue Liver

# Build allpairs table for coloc (set p-threshold to 0)
hvantk mktable eqtl \
  --raw-input /data/gtex_v11/allpairs/Liver/ \
  --output-ht eqtl_allpairs_liver.ht \
  --source gtex_v11 \
  --tissue Liver \
  --p-threshold 0
```

**GTEx v8** (TSV format):

```bash
# Download from GTEx Portal v8 archive
hvantk mktable eqtl \
  --raw-input /data/gtex_v8/Liver.v8.signif_variant_gene_pairs.txt.gz \
  --output-ht eqtl_liver_v8.ht \
  --source gtex_v8
```

**eQTLGen** (blood eQTLs):
URL: https://www.eqtlgen.org/cis-eqtls.html

```bash
# Download cis-eQTL full results (~2 GB)
hvantk mktable eqtl \
  --raw-input /data/eqtlgen/cis-eQTLs_full.txt.gz \
  --output-ht eqtl_blood.ht \
  --source eqtlgen
```

### Fang et al. (2025) pQTL data

Protein quantitative trait loci from Fang et al. (2025), covering 5 tissues (Colon, Heart, Liver, Lung, Thyroid). Space-delimited allpairs format with columns: `gene_name SNP CHR BP A1 NMISS BETA STAT P`. SE is derived as `|BETA/STAT|` (rows with `STAT = 0` are filtered out).

URL: Contact authors or GTEx Portal supplementary data.

> **Note:** Fang pQTL data uses gene symbols. Provide a `--gene-map-ht` (Ensembl gene table) for symbol → Ensembl ID mapping.

```bash
# Build pQTL table with gene mapping
hvantk mktable pqtl \
  --raw-input /data/fang_pqtl/Liver_allpairs.txt.gz \
  --output-ht pqtl_liver.ht \
  --source gtex_fang \
  --tissue Liver \
  --gene-map-ht ensembl_gene.ht \
  --p-threshold 5e-8

# Allpairs for coloc (omit p-threshold)
hvantk mktable pqtl \
  --raw-input /data/fang_pqtl/Liver_allpairs.txt.gz \
  --output-ht pqtl_allpairs_liver.ht \
  --source gtex_fang \
  --tissue Liver \
  --gene-map-ht ensembl_gene.ht
```

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
