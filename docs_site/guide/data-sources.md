# Data Sources

This page covers all annotation and expression data sources supported by hvantk: what they are, where to get them, and how to build datasets from the raw data. Sources are split into two categories: those with **built-in downloaders** (automated) and those that require **manual download** (too large, license-gated, or fragile URLs).

Every build on this page goes through the unified `hvantk reprocess <plugin>:<dataset>` entry point — see the [Usage Guide](usage.md#1-build-a-dataset-with-hvantk-reprocess) for the orchestration pattern, lifecycle flags, and `--plugin-arg KEY=VALUE` conventions used below.

## File format note

Downloaded `.gz` files may be standard gzip (single-threaded in Hail) rather than BGZF (parallel). Pre-convert before building:

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
| gnomAD constraint | `hvantk download gnomad-metrics` | ~4.6 MB (v2.1.1) / ~82 MB (v4.0) |
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
# Single command: download into data/hgnc/ then build the table
hvantk reprocess hgnc:lookup \
  --raw-dir data/hgnc/ \
  --output hgnc.ht

# Or, if you already downloaded the TSV into data/hgnc/:
hvantk reprocess hgnc:lookup \
  --raw-dir data/hgnc/ \
  --output hgnc.ht \
  --skip-download
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

These sources are too large, require license acceptance, or have complex download procedures. Follow the instructions below, place the raw file(s) in a per-source directory, then build via `hvantk reprocess <plugin>:<dataset> --skip-download`.

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
# Place the concatenated BGZF file in data/dbnsfp/ (so the builder sees it),
# then build. dbNSFP has no plugin downloader; --skip-download is required.
hvantk reprocess dbnsfp:variants \
  --raw-dir data/dbnsfp/ \
  --output dbnsfp.ht \
  --skip-download
```

> **Note:** dbNSFP's builder reads BGZF input. If you only have a single combined `.gz`, pre-convert with `hvantk utils convert-bgz dbNSFP4.9a_variant.gz` before building.

### gnomAD constraint metrics

Per-gene constraint metrics (pLI, oe_lof / LOEUF, missense Z) from gnomAD. The
tables are small and public, so hvantk ships a downloader
(`hvantk download gnomad-metrics`; it is also in the built-in-downloader list above).

- **v2.1.1** (GRCh37, ~4.6 MB) — the default and the table hvantk standardises on
  (keyed by `gene_id`).
- **v4.0** (GRCh38, ~82 MB) — newer; per-transcript rows, dotted column names, and
  **no `gene_id`** column, so build with `--plugin-arg key=transcript`. gnomAD did
  not re-release constraint for v4.1.

**Download + build** (end-to-end, defaults to v2.1.1 by-gene):

```bash
hvantk reprocess gnomad-metrics:metrics \
  --raw-dir data/gnomad_metrics/ \
  --output gnomad_metrics.ht

# v4.0 (GRCh38) instead:
hvantk download gnomad-metrics --version v4.0 \
  --output data/gnomad_metrics/gnomad.v4.0.constraint_metrics.tsv
hvantk reprocess gnomad-metrics:metrics --skip-download \
  --raw-dir data/gnomad_metrics/ --output gnomad_metrics_v4.ht \
  --plugin-arg key=transcript
```

### INSIDER interactome (~100 MB)

Protein-protein interaction sites from the INSIDER database.
URL: http://interactomeinsider.yulab.org/downloads.html

**Download**: Visit http://interactomeinsider.yulab.org/downloads.html and download the interaction site BED file.

**Build**:

```bash
# Place insider_interaction_sites.bed.bgz in data/insider/ then:
hvantk reprocess insider:variants \
  --raw-dir data/insider/ \
  --output interactome.ht \
  --skip-download
```

### Ensembl gene annotations (~60 MB GTF)

The canonical per-gene table (`ensembl-gene:structure`): gene ID, gene name, biotype,
coordinates, CDS length, coding-exon count, transcript count, and MANE Select — parsed from
the pinned-release Ensembl GTF. The release is pinned in
`hvantk/resources/ensembl_release.py`; the same pin governs the PTM coordinate mapper.

**Download**: the release-pinned GTF from the Ensembl FTP
(https://ftp.ensembl.org/pub/), or via the plugin's downloader.

**Build**:

```bash
hvantk reprocess ensembl-gene:structure \
  --raw-dir data/ensembl/ \
  --output ensembl_structure.ht \
  --skip-download
```

### GeVIR (~1-2 MB)

Gene-level intolerance-to-variation ranks (GeVIR and VIRLoF percentiles) for
19,361 protein-coding genes, keyed by Ensembl gene_id. Not a variant-level
pathogenicity score. Abramovs, Brass & Tassabehji, 2020, Nature Genetics
52(1):35-39 (PMID 31873297, DOI 10.1038/s41588-019-0560-2).
URL: https://www.nature.com/articles/s41588-019-0560-2

**Download**: Small supplementary table from the Nature Genetics publication
(https://www.nature.com/articles/s41588-019-0560-2) or the authors' repository
(https://github.com/gevirank/gevir). At ~1-2 MB with a stable, public URL, GeVIR
qualifies for a real downloader under the framework in CLAUDE.md — a recommended
follow-up (not yet implemented).

**Build**:

```bash
# Place gevir_metrics.tsv.bgz in data/gevir/ then:
hvantk reprocess gevir:metrics \
  --raw-dir data/gevir/ \
  --output gevir.ht \
  --skip-download
```

### COSMIC Cancer Gene Census

Gene-level cancer annotations from the COSMIC Cancer Gene Census.
URL: https://cancer.sanger.ac.uk/census

**Download**: Requires COSMIC account. Download the Cancer Gene Census TSV from the COSMIC website.

**Build**:

```bash
# Place cancer_gene_census.tsv in data/cosmic_cgc/ then:
hvantk reprocess cosmic-cgc:submissions \
  --raw-dir data/cosmic_cgc/ \
  --output cosmic_cgc.ht \
  --skip-download
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
# Navigate to GTEx Portal → Downloads → Adult GTEx → QTL → eQTL → Significant pairs.
# Place per-tissue files under /data/gtex_v11/signif_pairs/ then:

# Build significant-pairs table
hvantk reprocess gtex-eqtl:eqtls \
  --raw-dir /data/gtex_v11/signif_pairs/ \
  --output eqtl_liver.ht \
  --skip-download \
  --plugin-arg source=gtex_v11 \
  --plugin-arg tissue=Liver

# Build allpairs table for coloc (set p-threshold to 0)
hvantk reprocess gtex-eqtl:eqtls \
  --raw-dir /data/gtex_v11/allpairs/Liver/ \
  --output eqtl_allpairs_liver.ht \
  --skip-download \
  --plugin-arg source=gtex_v11 \
  --plugin-arg tissue=Liver \
  --plugin-arg p_threshold=0
```

**GTEx v8** (TSV format):

```bash
# Place Liver.v8.signif_variant_gene_pairs.txt.gz under /data/gtex_v8/ then:
hvantk reprocess gtex-eqtl:eqtls \
  --raw-dir /data/gtex_v8/ \
  --output eqtl_liver_v8.ht \
  --skip-download \
  --plugin-arg source=gtex_v8
```

**eQTLGen** (blood eQTLs):
URL: https://www.eqtlgen.org/cis-eqtls.html

```bash
# Place cis-eQTLs_full.txt.gz under /data/eqtlgen/ then:
hvantk reprocess gtex-eqtl:eqtls \
  --raw-dir /data/eqtlgen/ \
  --output eqtl_blood.ht \
  --skip-download \
  --plugin-arg source=eqtlgen
```

### Fang et al. (2025) pQTL data

Protein quantitative trait loci from Fang et al. (2025), covering 5 tissues (Colon, Heart, Liver, Lung, Thyroid). Space-delimited allpairs format with columns: `gene_name SNP CHR BP A1 NMISS BETA STAT P`. SE is derived as `|BETA/STAT|` (rows with `STAT = 0` are filtered out).

URL: Contact authors or GTEx Portal supplementary data.

> **Note:** Fang pQTL data uses gene symbols. Pass an Ensembl gene-table path via `--plugin-arg hgnc_ht=<path>` for symbol → Ensembl ID mapping (the builder uses an HGNC-style lookup table).

```bash
# Build pQTL table with gene mapping
# Place Liver_allpairs.txt.gz under /data/fang_pqtl/ then:
hvantk reprocess pqtl:metrics \
  --raw-dir /data/fang_pqtl/ \
  --output pqtl_liver.ht \
  --skip-download \
  --plugin-arg source=gtex_fang \
  --plugin-arg tissue=Liver \
  --plugin-arg hgnc_ht=ensembl_gene.ht \
  --plugin-arg p_threshold=5e-8

# Allpairs for coloc (omit p_threshold to keep all variants)
hvantk reprocess pqtl:metrics \
  --raw-dir /data/fang_pqtl/ \
  --output pqtl_allpairs_liver.ht \
  --skip-download \
  --plugin-arg source=gtex_fang \
  --plugin-arg tissue=Liver \
  --plugin-arg hgnc_ht=ensembl_gene.ht
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
