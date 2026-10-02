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
| Ensembl gene annotations | `hvantk download ensembl-structure` | ~60 MB |
| UniProt PTM sites | `hvantk download uniprot-ptm` | varies (REST API query) |
| CPTAC phosphoproteomics | `hvantk download cptac-phospho` | varies by cancer type |
| PeptideAtlas phospho | `hvantk download peptideatlas-phospho` | ~549 MB |
| 1000 Genomes sample panel | *(no standalone command; see [1000 Genomes](#1000-genomes-nygcccdg-high-coverage-callset))* | ~55 KB |

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
hvantk download hgnc --output data/hgnc/hgnc_complete_set.txt
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

Most of these sources need manual acquisition — too large, license-gated, or a fragile URL
— and for those, follow the instructions below, place the raw file(s) in a per-source
directory, then build via `hvantk reprocess <plugin>:<dataset> --skip-download`. A few
members of this section (gnomAD constraint, Ensembl gene annotations, UniProt PTM sites)
ship a `hvantk download` command instead and are grouped here only because they carry a
license, size, or versioning caveat worth reading before relying on them.

### dbNSFP (~50 GB)

A database of functional prediction scores for human missense variants.
URL: https://www.dbnsfp.org/

**Download**: Academic access is a two-step form process: (1) register with an institutional email through a Google Form and receive an access code once the domain is verified; (2) request the download links with that same email and code through a second form. The academic release is distributed under the CC BY-NC-ND 4.0 licence, as a ~50 GB ZIP of per-chromosome variant tables plus the gene table and the `search_dbNSFP` program; it is also searchable online through an "Academic Portal". Commercial licensing goes through the maintainers. Start at the releases page, https://www.dbnsfp.org/releases/, or the download page, https://www.dbnsfp.org/download/. (Verified 2026-09-28: the legacy landing page, https://sites.google.com/site/jpopgen/dbNSFP, is frozen at v4.9 (2024-08-08) and now tells readers to access the current dbNSFP at dbNSFP.org; its S3 and Box archive links are dead.)

**Pre-processing**: dbNSFP is distributed as per-chromosome `.gz` files (standard gzip, not BGZF). The builder expects a **single combined file**, so concatenate and BGZF-compress first:

```bash
# Concatenate per-chromosome files into a single BGZF file
# (header is taken from chr1; remaining files skip the header line)
head -1 <(zcat dbNSFP<version>_variant.chr1.gz) > /tmp/dbnsfp_header.txt
(cat /tmp/dbnsfp_header.txt && for f in dbNSFP<version>_variant.chr*.gz; do zcat "$f" | tail -n +2; done) \
  | bgzip -@ 4 > dbnsfp_variant.bgz
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

> **Note:** dbNSFP's builder reads BGZF input. If you only have a single combined `.gz`, pre-convert with `hvantk utils convert-bgz dbNSFP<version>_variant.gz` before building.

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

### INSIDER interactome (~1.2 GB genomic BED; ~49 MB pair table)

Protein-protein interaction interface residues from the INSIDER database.
URL: http://interactomeinsider.yulab.org/downloads.html

INSIDER ships **two** products, and hvantk builds a dataset from each:

| dataset | file | size | direct URL |
| --- | --- | --- | --- |
| `insider:variants` | `Whole_Human_Interactome_Interface_hg38.bed` | ~1.17 GB | `http://interactomeinsider.yulab.org/bed/all.bed` |
| `insider:interfaces` | `H_sapiens_interfacesALL.txt` | ~49 MB | `http://interactomeinsider.yulab.org/downloads/interfacesALL/H_sapiens_interfacesALL.txt` |

**Download**: the downloads page carries no links in its markup, so use the direct
URLs above (they are also recorded in the plugin catalog — `hvantk catalog show
INSIDER_v1.0`). Both are served over plain HTTP; the site has no HTTPS listener.
The BED is >1 GB, so acquisition is manual per the downloader framework in CLAUDE.md.

**Build**:

```bash
# Place insider_interaction_sites.bed.bgz in data/insider/ then:
hvantk reprocess insider:variants \
  --raw-dir data/insider/ \
  --output interactome.ht \
  --skip-download
```

INSIDER also ships a **gene-level** dataset, `insider:interfaces`, built from
`H_sapiens_interfacesALL.txt` (~49 MB) rather than the BED file. Where `insider:variants` is
locus-keyed (interaction sites as genomic intervals), `insider:interfaces` is keyed by
`uniprot_id` and summarises each protein's interactome: `n_partners`,
`n_partners_experimental`, `n_partners_predicted` and `n_interface_residues`. The two are
complementary, not alternatives — use `interfaces` when the question is "how connected is
this gene", and `variants` when it is "does this variant fall in an interaction site". Unlike
`insider:variants`, this file is well under 1 GB, but `insider:interfaces` declares no
`acquisition` block and no `lifecycle.download` either — no downloader exists yet for it
(tracked in #386), so `--skip-download` is required below for the same reason as `gevir` and
`gwas-catalog`.

```bash
# Place H_sapiens_interfacesALL.txt in data/insider/ then:
hvantk reprocess insider:interfaces \
  --raw-dir data/insider/ \
  --output insider_interfaces.ht \
  --skip-download
```

### Ensembl gene annotations (~60 MB GTF)

The canonical per-gene table (`ensembl-gene:structure`): gene ID, gene name, biotype,
coordinates, CDS length, coding-exon count, transcript count, and MANE Select — parsed from
the pinned-release Ensembl GTF. The release is pinned in
`hvantk/resources/ensembl_release.py`; the same pin governs the PTM coordinate mapper. The
file is small and public, so hvantk ships a downloader (`hvantk download ensembl-structure`;
it is also in the built-in-downloader list above) that fetches exactly that pinned release
rather than whatever the FTP's "current" alias points at.

**Download + build** (end-to-end, fetches the pinned release automatically):

```bash
hvantk reprocess ensembl-gene:structure \
  --raw-dir data/ensembl/ \
  --output ensembl_structure.ht
```

**Manual download**: the same release-pinned GTF from the Ensembl FTP
(https://ftp.ensembl.org/pub/), placed in `data/ensembl/`, then built with `--skip-download`:

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

**Download**: the metric table is **Supplementary Table 2** of the Nature Genetics
paper, served as the article's MOESM3 object:

```
https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0560-2/MediaObjects/41588_2019_560_MOESM3_ESM.xlsx
```

That is a ~10.3 MB `.xlsx` workbook (only MOESM3 of the six supplementary slots is
public; the rest return 403). The builder reads a bgzipped TSV, so extract sheet
`table_2` and BGZF-compress it before building.

> **Note:** the authors' repository at https://github.com/gevirank/gevir ships the
> **analysis code only** — its `tables/` directory holds a placeholder file — so it
> is not a source for the metric table. Earlier revisions of this guide pointed
> there.

A real downloader would have to do the extract-and-convert step, not just fetch the
URL, so it is more than the usual thin wrapper — tracked as a follow-up in #386.

**Build**:

```bash
# Place gevir_metrics.tsv.bgz in data/gevir/ then:
hvantk reprocess gevir:metrics \
  --raw-dir data/gevir/ \
  --output gevir.ht \
  --skip-download
```

### GWAS Catalog (~69.5 MB zip, ~601 MB unzipped)

The NHGRI-EBI GWAS Catalog full-associations TSV — SNP-trait associations with
study provenance, p-value, effect size, and risk-allele frequency, built into a
Hail Table keyed by `(locus, alleles)` with a sentinel ALT. v1.0 schema (34 columns).
URL: https://www.ebi.ac.uk/gwas/

**Download**: a periodic full dump, zipped, from EBI's FTP, released every 1-3 weeks:

```
https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip
```

Unzip to get `gwas-catalog-download-associations-v1.0-full.tsv` (verified live
2026-09-30: zip 69,506,531 bytes, TSV 600,996,055 bytes). No `lifecycle.download` is
declared yet for this plugin — see #386.

**Build**:

```bash
# Place the unzipped TSV in data/gwas_catalog/ then:
hvantk reprocess gwas-catalog:associations \
  --raw-dir data/gwas_catalog/ \
  --intermediate data/gwas_catalog/gwas-catalog-download-associations-v1.0-full.tsv \
  --skip-parse \
  --skip-download \
  --output gwas_catalog.ht
```

> **Note:** `--intermediate <file> --skip-parse` is required, not optional, even though
> this dataset declares no `lifecycle.parse` — the builder needs the exact TSV path, and
> `--raw-dir` can only be a directory. `--skip-download` is also required — no
> `lifecycle.download` is declared. See `hvantk/skills/gwas_catalog/SKILL.md` §6.

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

### AlphaGenome variant effect predictions

Per-variant deep-learning effect predictions (expression, chromatin, and other molecular
phenotypes) from Google DeepMind's AlphaGenome model (`alphagenome:predictions`), fetched
live through a credentialed API rather than read from any downloadable file. SDK:
https://pypi.org/project/alphagenome/

**Access**: provision API credentials yourself — `api.key` in a config YAML, or the
`ALPHAGENOME_API_KEY` environment variable; hvantk does not provision them. The config YAML
needs `api` and `ontology` sections (ontology terms plus output types such as `RNA_SEQ`,
`CHROMATIN`); see `hvantk/skills/alphagenome/tests/testdata/alphagenome_config.yaml` for the
shape.

**Build**:

```bash
# Place a variant Hail Table (or a directory containing one) at data/alphagenome/variants.ht
hvantk reprocess alphagenome:predictions \
  --raw-dir data/alphagenome/variants.ht \
  --output alphagenome.ht \
  --skip-download \
  --plugin-arg config_path=path/to/alphagenome_config.yaml
```

> **Note:** `--skip-download` is not strictly required here — `acquisition.mode: byo`
> already makes skipping implicit — but passing it stays legal and matches the rest of
> this section. `--raw-dir` must be a directory (a Hail Table satisfies this; a plain
> `.tsv` does not), since the plugin declares no `lifecycle.parse` and `reprocess` forwards
> `--raw-dir` straight to the builder. Each build issues live, billed API calls with no
> cross-run resume, so budget accordingly and prefer chunking a large variant set.

### MSigDB gene sets

Molecular Signatures Database gene sets (`msigdb:genesets`), built from a collection's GMT
file (e.g. C2 canonical pathways). The format is the same across collections, so the
builder works for any MSigDB GMT.
URL: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp

**Download**: requires a free account and a login + license click-through, so no
downloader is shipped. Download and unzip the collection's GMT.

**Build**:

```bash
# Place the unzipped .gmt in data/msigdb/ then:
hvantk reprocess msigdb:genesets \
  --raw-dir data/msigdb/ \
  --output msigdb.ht \
  --skip-download
```

### UniProt PTM Sites

Curated post-translational modification sites (phosphorylation, ubiquitination, acetylation, etc.) for reviewed human proteins from UniProt/Swiss-Prot.
URL: https://www.uniprot.org/

This is a live REST query rather than a large or license-gated file, so hvantk ships a
downloader (`hvantk download uniprot-ptm`; it is also in the built-in-downloader list
above). The `uniprot-ptm:sites` plugin dataset builds from a *genome-mapped* TSV, not the
raw UniProt download, so `hvantk reprocess uniprot-ptm:sites` alone cannot complete the
coordinate-mapping step — use `hvantk ptm build` below, which downloads the raw TSV, maps
it against the Ensembl GTF, and builds the Hail Table in one command. For pre-download or
manual acquisition:

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

### CPTAC proteomics

Clinical Proteomic Tumor Analysis Consortium (CPTAC) tumor proteomics, accessed via the
[`cptac` Python package](https://pypi.org/project/cptac/) (not direct HTTP) rather than a
single downloadable file.
URL: https://proteomics.cancer.gov/programs/cptac

**`cptac:phospho`** (phosphoproteomics, one AnnData per cancer type) ships a downloader
(`hvantk download cptac-phospho`; it is also in the built-in-downloader list above):

```bash
hvantk reprocess cptac:phospho \
  --raw-dir data/cptac/ \
  --output cptac_phospho_brca.ht \
  --plugin-arg cancer_type=brca
```

**`cptac:expression`** (protein expression) has no downloader and `hvantk reprocess
cptac:expression` is not wired at all: the manifest declares no `lifecycle.parse` to turn
a raw download into the `{"expression": <tsv>, "metadata": <tsv>}` mapping the builder
needs. Stage the expression and metadata TSVs yourself and call the builder directly —
see `hvantk/skills/cptac/expression/SKILL.md` §6-7 for the full workflow
(`build_cptac_expression(parsed_input, ctx, gene_id_col=..., sample_id_col=...,
expression_col=...)`); an automated downloader is tracked in #386.

### PeptideAtlas phospho sites

Human phosphorylation sites from PeptideAtlas (`peptideatlas:phospho`), built from the
PeptideAtlas relational TSV dump.
URL: https://peptideatlas.org/

```bash
hvantk reprocess peptideatlas:phospho \
  --raw-dir data/peptideatlas/ \
  --output peptideatlas_phospho.ht
```

The upstream zip (e.g. `atlas_build_606.tsv.zip`) is ~549 MB. `hvantk download
peptideatlas-phospho` (also in the built-in-downloader list above) fetches it, and the
plugin's own `lifecycle.parse` stage extracts the intermediate TSV automatically — no
`--skip-parse` or `--intermediate` needed.

### 1000 Genomes (NYGC/CCDG high-coverage callset)

Per-chromosome high-coverage VCFs (NYGC/CCDG) plus IGSR sample metadata (population,
super-population, sex).

**`onek-genomes:variants`** (the VCFs) is `byo` (reason: size — multi-GB per chromosome),
so no downloader is shipped. Mirror:
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/

```bash
# Place per-chromosome *_chr<N>.filtered.shapeit2-duohmm-phased.vcf.gz (+ .tbi) in
# data/onek_genomes/ then:
hvantk reprocess onek-genomes:variants \
  --raw-dir data/onek_genomes/ \
  --output onek_genomes_variants.ht \
  --skip-download
```

**`onek-genomes:samples`** (the IGSR sample panel, ~55 KB) declares a `lifecycle.download`,
but the plugin has no `cli.py` and no `cli:` block, so no standalone `hvantk download`
command exists for it — the download only runs as part of `reprocess`:

```bash
hvantk reprocess onek-genomes:samples \
  --raw-dir data/onek_genomes/ \
  --output onek_genomes_samples.ht
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
