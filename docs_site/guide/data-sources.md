# Data Sources

This page covers all annotation and expression data sources supported by hvantk: what they are, where to get them, and how to build datasets from the raw data. Sources are split into two categories: those with **built-in downloaders** (automated) and those that require **manual download** (too large, license-gated, or fragile URLs).

Every build on this page goes through the unified `hvantk reprocess <plugin>:<dataset>` entry point — see [Building datasets](usage.md#build-a-dataset) for the orchestration pattern, lifecycle flags, and `--plugin-arg KEY=VALUE` conventions used below.

## Dataset summary

Every dataset `hvantk` can build, grouped by how you get the raw data. `Backend` is the
builder's output format (`hail` writes a Hail Table, `.ht`; `anndata` writes an AnnData,
`.h5ad`; `pandas` writes a Parquet file, `.parquet`). Licence/terms notes are what the
plugin's own catalog entry, its `SKILL.md` or
[THIRD_PARTY_DATA.md](https://github.com/bigbio/hvantk/blob/main/THIRD_PARTY_DATA.md)
states; the last quotes the source of each licence it lists. "Not stated in plugin docs"
means the terms are not documented in this repository, not that the data is
unrestricted. Check the source's own terms before you redistribute anything you build.

### Built-in downloaders (16 datasets)

| Dataset | Type | Backend | Licence / terms |
|---|---|---|---|
| `clinvar:variants` | AnnotationTable | hail | Public domain |
| `clingen:gene-disease` | AnnotationTable | hail | CC BY 4.0 |
| `gencc:submissions` | AnnotationTable | hail | not stated in plugin docs |
| `hgnc:lookup` | AnnotationTable | hail | CC0 |
| `gnomad-metrics:metrics` | AnnotationTable | hail | CC0 |
| `ucsc-cellbrowser:default` | ExpressionMatrix | anndata | varies per collection; not centrally stated |
| `expression-atlas:dataset` | ExpressionMatrix | anndata | varies per accession (CC BY 4.0 for several; see `hvantk catalog show <accession>`) |
| `ensembl-gene:structure` | AnnotationTable | hail | No restrictions (Ensembl data; third-party constraints may apply) |
| `uniprot-ptm:sites` | AnnotationTable | hail | not stated in plugin docs |
| `cptac:expression` | ExpressionMatrix | anndata | not stated in plugin docs (via the `cptac` PyPI package) |
| `cptac:phospho` | ExpressionMatrix | anndata | not stated in plugin docs (via the `cptac` PyPI package) |
| `peptideatlas:phospho` | AnnotationTable | pandas | not stated in plugin docs |
| `insider:interfaces` | AnnotationTable | hail | Not stated (checked 2026-10-13; see `SKILL.md`) |
| `gwas-catalog:associations` | AnnotationTable | hail | EBI Services Terms of Use (CC0 covers only the separate summary statistics) |
| `gevir:metrics` | AnnotationTable | hail | Academic use |
| `onek-genomes:samples` | AnnotationTable | hail | No restriction; acknowledgement required (downloads only through `reprocess`; there is no `hvantk download` command for it) |

### Bring your own (10 datasets)

| Dataset | Type | Backend | Reason | Licence / terms |
|---|---|---|---|---|
| `alphagenome:predictions` | AnnotationTable | hail | credentialed | AlphaGenome Output Terms of Use: non-commercial |
| `cosmic-cgc:submissions` | AnnotationTable | hail | credentialed | COSMIC account required; licence not stated in plugin docs |
| `dbnsfp:variants` | AnnotationTable | hail | size | CC BY-NC-ND 4.0 (academic); commercial licensing via the maintainers |
| `gtex-eqtl:eqtls` | AnnotationTable | hail | unstable-url | Open access, with GTEx's publication embargo and citation rules |
| `insider:variants` | AnnotationTable | hail | size | Not stated (checked 2026-10-13; see `SKILL.md`) |
| `msigdb:genesets` | AnnotationTable | hail | credentialed | CC Attribution 4.0 (academic) / commercial licence required; BioCarta and legacy KEGG sets are excluded from the CC BY grant |
| `onek-genomes:variants` | VariantMatrix | hail | size | No restriction; acknowledgement required |
| `pqtl:metrics` | AnnotationTable | hail | publication-only | Preprint supplementary material (Fang et al. 2025); no licence stated |
| `ucsc-cellbrowser:adult-ctx` | ExpressionMatrix | anndata | size | varies (parent collection); not centrally stated |
| `ucsc-cellbrowser:dev-ctx` | ExpressionMatrix | anndata | size | varies (parent collection); not centrally stated |

## Downloader or documented acquisition

A dataset's manifest declares one of three acquisition states — per *dataset*, not
per provider, since a multi-dataset plugin can mix them (`onek-genomes` ships one BYO
dataset and one auto-downloaded one):

- **`lifecycle.download` is declared** — a real downloader exists; `hvantk download
  <source>` or `hvantk reprocess <plugin>:<dataset>` fetches the raw file.
- **`acquisition: {mode: byo, reason: ...}` is declared** — acquisition cannot be
  automated, so it is documented instead (`reason` is one of `size`, `credentialed`,
  `license`, `publication-only`, `unstable-url`). The user supplies `--raw-dir`
  themselves; `byo` and `lifecycle.download` are mutually exclusive — a manifest
  cannot claim both that fetching is impossible and ship a function that fetches.
  `reprocess` skips the download stage automatically for a BYO dataset — passing
  `--skip-download` yourself is legal but not required.
- **Neither block is declared** — automatable in principle, but no downloader has
  been written yet; `hvantk reprocess` needs `--skip-download` until one lands.

(An explicit `acquisition.mode: download` with no `lifecycle.download` is a schema
error: it asserts a downloader exists without shipping one.)

The decision between the first two states: build a downloader when the raw file is
**under 1 GB**, **publicly accessible** (no licence gate), **high-value**, and served
from a **stable, predictable URL**. Document acquisition instead when the file is
**over 1 GB**, needs **licence acceptance or manual steps**, or the URL structure is
**complex or fragile** (FTP directory listings, versioned mirrors). The third state
means a downloader could be written in principle but has not been yet; nothing on
this page is currently in that state.

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
| Ensembl gene annotations | `hvantk download ensembl-structure` | ~64 MB |
| UniProt PTM sites | `hvantk download uniprot-ptm` | varies (REST API query) |
| CPTAC protein expression | `hvantk download cptac-expression` | varies by cancer type |
| CPTAC phosphoproteomics | `hvantk download cptac-phospho` | varies by cancer type |
| PeptideAtlas phospho | `hvantk download peptideatlas-phospho` | ~549 MB |
| INSIDER interfaces | `hvantk download insider-interfaces` | ~49 MB |
| GWAS Catalog | `hvantk download gwas-catalog` | ~70 MB zip (~601 MB unzipped TSV) |
| GeVIR gene metrics | `hvantk download gevir` | ~10.3 MB xlsx (~0.6 MB derived TSV) |
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

**Build**:

```bash
# Download + build in one step
hvantk reprocess clinvar:variants \
  --raw-dir data/clinvar/ \
  --output clinvar.ht

# From a VCF already in --raw-dir
hvantk reprocess clinvar:variants \
  --raw-dir data/clinvar/ \
  --output clinvar.ht \
  --skip-download \
  --plugin-arg reference_genome=GRCh38
```

### ClinGen

```bash
# Download today's ClinGen Gene-Disease Validity snapshot
# Output: Clingen-Gene-Disease-Summary-<YYYY-MM-DD>.csv
hvantk download clingen --output-dir data/clingen

# Check download availability
hvantk download clingen --list-versions
```

**Build**:

```bash
hvantk reprocess clingen:gene-disease \
  --raw-dir data/clingen/ \
  --output clingen.ht \
  --plugin-arg min_classification=Strong
```

### GenCC

GenCC (Gene Curation Coalition) aggregates gene-disease validity assertions from 12+ submitting organizations (ClinGen, PanelApp, G2P, Orphanet, etc.).

```bash
# Download today's GenCC submissions snapshot
hvantk download gencc --output-dir data/gencc

# Check download availability
hvantk download gencc --list-versions
```

**Build**:

```bash
hvantk reprocess gencc:submissions \
  --raw-dir data/gencc/ \
  --output gencc.ht \
  --plugin-arg min_classification=Strong
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

The UCSC Cell Browser catalog lists 270 datasets. About half are **collections** (groups
of related datasets with no expression matrix at the top level). Use `--list-datasets`
and `--search` to discover downloadable datasets.

```bash
# Discover available datasets
hvantk download ucsc --list-datasets

# Search by name, organism, or tissue (expands collections to show children)
hvantk download ucsc --list-datasets --search heart
hvantk download ucsc --list-datasets --search pancreas

# Download a leaf dataset directly
hvantk download ucsc --dataset adultPancreas --output-dir data/ucsc

# Download a child dataset from a collection (use the full path)
hvantk download ucsc --dataset hoc/all-heart --output-dir data/ucsc
```

> **Note:** Collection names (e.g., `hoc`) cannot be downloaded directly — they
> contain no expression matrix. Use `--search` to find child dataset paths like
> `hoc/all-heart`, then download those.

**Build** (`ucsc-cellbrowser:default`, one named collection per build):

```bash
# Download + build in one step (dataset forwards to the download stage)
hvantk reprocess ucsc-cellbrowser:default \
  --raw-dir data/ucsc/adultPancreas \
  --output adult_pancreas.h5ad \
  --plugin-arg dataset=adultPancreas

# From files already downloaded with `hvantk download ucsc` above
hvantk reprocess ucsc-cellbrowser:default \
  --raw-dir data/ucsc/adultPancreas \
  --output adult_pancreas.h5ad \
  --skip-download
```

**`adult-ctx` and `dev-ctx`** are per-class summaries of the `adult-ctx-meta-atlas` and
`dev-ctx-meta-atlas` parent collections, not their own UCSC collections, so neither
declares `lifecycle.download` (`acquisition.mode: byo`, reason `size` — the parent's raw
per-cell matrix is multi-GB). Each needs an `expression_matrix.tsv` (a `gene` column plus
one column per class) and a `metadata.tsv` (one row per class, with `n_cells`) in
`--raw-dir`; `hvantk expression summarize-ucsc` builds that summary as an `.h5ad` from the
parent collection, so exporting it back to this TSV pair takes a manual step:

```bash
hvantk reprocess ucsc-cellbrowser:adult-ctx \
  --raw-dir data/ucsc/adult_ctx/ \
  --output adult_ctx_summary.h5ad \
  --skip-download

hvantk reprocess ucsc-cellbrowser:dev-ctx \
  --raw-dir data/ucsc/dev_ctx/ \
  --output dev_ctx_summary.h5ad \
  --skip-download
```

### Expression Atlas

```bash
# Download one accession's bulk RNA-seq experiment
hvantk download expression-atlas --accession E-MTAB-6798 --download-path data/expression_atlas

# Download + build in one call (reprocess forwards --plugin-arg accession=... to the
# download stage; the build stage ignores it and locates the two input files within
# --raw-dir)
hvantk reprocess expression-atlas:dataset --raw-dir data/expression_atlas \
  --output E-MTAB-6798.h5ad --plugin-arg accession=E-MTAB-6798
```

`reprocess` scopes the expression matrix to the accession named by the raw dir's SDRF
file and raises if no TPM file matches it (an older raw dir's leftover file, or two
downloads merged by hand, must not silently combine mismatched expression + metadata).
The unfiltered download above fetches every file in the experiment's FTP directory, so
the gene-level and the transcript-level TPM file can both be present. `reprocess` builds
from the one whose columns are the SDRF's runs, or from the **gene-level** one when both
qualify. For E-MTAB-6798 that is the **transcript-level** file (`var` keyed by the
transcript id): its gene-level file has one column per assay group (`g1` … `g94`), not
per run. It raises when no file has the runs as columns.

A transcript-level build needs `--plugin-arg level=gene` to sum its values per gene and
key `var` by the gene id instead — before the result can be used with
`hvantk expression markers`/`summarize` or the PTM constraint expression loader, which
read `var_names` as the gene id and refuse a transcript-keyed matrix. A gene-level build
already works there without it. `hvantk annotate prepare`'s
matrix source type refuses a transcript-keyed matrix too, but needs more than
`level=gene` either way: it
resolves `var_names` as HGNC gene symbols, not Ensembl gene ids, so an Expression Atlas
build (gene- or transcript-level) needs a symbol-keyed summary before it is usable
there.

## Manual download sources

Most of these sources need manual acquisition — too large, license-gated, or a fragile URL
— and for those, follow the instructions below, place the raw file(s) in a per-source
directory, then build via `hvantk reprocess <plugin>:<dataset>` (a BYO dataset skips the
download stage automatically; `--skip-download` is legal but not required). Some
members of this section ship a downloader instead (GeVIR, gnomAD constraint, Ensembl gene
annotations, UniProt PTM sites, CPTAC phospho and expression, PeptideAtlas phospho,
INSIDER interfaces, the GWAS Catalog and the 1000 Genomes sample metadata); they are grouped here because
they carry a license, size, versioning or build-path caveat worth reading before relying
on them.

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
# Place the concatenated BGZF file in data/dbnsfp/ (so the builder sees it), then
# build. dbNSFP declares acquisition.mode=byo (no plugin downloader exists), so
# reprocess skips the download stage automatically.
hvantk reprocess dbnsfp:variants \
  --raw-dir data/dbnsfp/ \
  --output dbnsfp.ht
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

**Download**: the downloads page carries no links in its markup, so both the plugin
catalog (`hvantk catalog show INSIDER_v1.0`) and the `insider:interfaces` downloader use
the direct URLs above. Both files are served over plain HTTP; the site has no HTTPS
listener. The BED is >1 GB, so `insider:variants` stays manual acquisition; the pair
table is small enough that `insider:interfaces` ships a real downloader instead:

```bash
hvantk download insider-interfaces --output data/insider/H_sapiens_interfacesALL.txt
```

**Build**:

```bash
# Place Whole_Human_Interactome_Interface_hg38.bed in data/insider/ then:
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
this gene", and `variants` when it is "does this variant fall in an interaction site".
Unlike `insider:variants`, this file is well under 1 GB, so it declares `lifecycle.download`
and `reprocess` fetches it automatically — no `--skip-download` needed:

```bash
hvantk reprocess insider:interfaces \
  --raw-dir data/insider/ \
  --output insider_interfaces.ht
```

### Ensembl gene annotations (~64 MB GTF)

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

### GeVIR (~10 MB workbook)

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
public; the rest return 403), not a thin URL fetch: the builder reads a bgzipped TSV,
so the downloader (`hvantk download gevir`; it is also in the built-in-downloader list
above) extracts sheet `table_2` itself (stdlib `zipfile` + `xml.etree`, no openpyxl
dependency) and BGZF-compresses it, including restoring the ~29 gene symbols Excel's
autocorrect silently turned into dates (e.g. `MARCH6`, `SEPT11`; see
`hvantk/skills/gevir/cli.py`'s module docstring for the full, externally-verified
restoration rule).

> **Note:** the authors' repository at https://github.com/gevirank/gevir ships the
> **analysis code only** — its `tables/` directory holds a placeholder file — so it
> is not a source for the metric table. Earlier revisions of this guide pointed
> there.

**Download + build** (end-to-end, fetches and derives the TSV automatically):

```bash
hvantk reprocess gevir:metrics \
  --raw-dir data/gevir/ \
  --output gevir.ht
```

**Manual download**: an already-derived `gevir_metrics_pmid31873297.tsv.bgz` placed in
`data/gevir/`, then built with `--skip-download`:

```bash
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

A periodic full dump, zipped, from EBI's FTP, released every 1-3 weeks. The zip is small
and public, so hvantk ships a downloader (`hvantk download gwas-catalog`; it is also in
the built-in-downloader list above) that fetches the current release from the FTP
`releases/latest/` alias and extracts the associations TSV it contains.

**Download + build** (end-to-end, fetches the current release automatically):

```bash
hvantk reprocess gwas-catalog:associations \
  --raw-dir data/gwas_catalog/ \
  --output gwas_catalog.ht
```

**Manual download**:

```
https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip
```

Unzip to get `gwas-catalog-download-associations-v1.0-full.tsv` (verified live
2026-09-30: zip 69,506,531 bytes, TSV 600,996,055 bytes), place it in `data/gwas_catalog/`,
then build with `--skip-download`:

```bash
hvantk reprocess gwas-catalog:associations \
  --raw-dir data/gwas_catalog/ \
  --output gwas_catalog.ht \
  --skip-download
```

### COSMIC Cancer Gene Census

Gene-level cancer annotations from the COSMIC Cancer Gene Census.
URL: https://cancer.sanger.ac.uk/census

**Download**: Requires COSMIC account. Download the Cancer Gene Census TSV from the COSMIC website.

**Build**:

```bash
# Place the downloaded export (e.g. cancer_gene_census.tsv.gz) in data/cosmic_cgc/ then:
hvantk reprocess cosmic-cgc:submissions \
  --raw-dir data/cosmic_cgc/ \
  --skip-download \
  --output cosmic_cgc.ht
```

### AlphaGenome variant-effect scores

Per-variant scores from Google DeepMind's AlphaGenome model (`alphagenome:predictions`) for
expression, splicing, chromatin and contact-map tracks. hvantk does not call AlphaGenome: it
ingests scores you have already produced with the AlphaGenome SDK, through the API or local
weights, and summarises them into one row per variant. SDK: https://pypi.org/project/alphagenome/
(documentation: https://www.alphagenomedocs.com).

**Licence**: AlphaGenome outputs fall under the AlphaGenome Output Terms of Use: non-commercial
use only, not for training machine-learning models, and any copy keeps a conspicuous notice
(see `hvantk/skills/alphagenome/tests/testdata/raw/alphagenome/NOTICE.md`).

**Produce the input**: score each variant with `score_variant(...)`, convert the result with
`variant_scorers.tidy_scores(...)`, cast `variant_id` and `scored_interval` to strings, and write
the DataFrame with `to_parquet()`. `hvantk/skills/alphagenome/SKILL.md` §2 has the snippet. Put
only these parquet files in one directory; every `*.parquet` file in it is read.

**Build**:

```bash
# Place the tidy_scores parquet files in data/alphagenome/ then:
hvantk reprocess alphagenome:predictions \
  --raw-dir data/alphagenome/ \
  --output alphagenome.ht
```

> **Note:** the table is keyed by `(locus, alleles)` and carries one struct per AlphaGenome
> scorer (the 19 of the SDK's recommended set): the track with the largest absolute quantile
> score, the largest absolute raw score and the number of rows. `--plugin-arg output_types=...`
> and `--plugin-arg ontology_curies=...` filter rows before they are summarised. A scorer the
> plugin does not know, or a parquet file without the tidy-scores columns, stops the build with
> an error. `--skip-download` is not needed: the dataset declares `acquisition.mode: byo`.

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

Curated post-translational modification sites for reviewed human proteins from
UniProt/Swiss-Prot. The downloader queries only UniProt's `ft_mod_res` ("Modified
residue") feature — phosphorylation, acetylation, methylation, hydroxylation and
similar single-residue marks. It does not cover modifications UniProt annotates under
a different feature type: no ubiquitination, SUMOylation or NEDDylation (`CROSSLNK`)
and no glycosylation (`CARBOHYD`).
URL: https://www.uniprot.org/

This is a live REST query rather than a large or license-gated file, so hvantk ships a
downloader (`hvantk download uniprot-ptm`; it is also in the built-in-downloader list
above). The `uniprot-ptm:sites` table is built from *genome-mapped* sites, not from the
raw UniProt download, and `hvantk reprocess uniprot-ptm:sites` runs every step: it
downloads the raw TSV, fetches the pinned Ensembl GTF into the same `--raw-dir` (kept
there, so later runs reuse it), maps the sites against it, and builds the Hail Table.
`hvantk ptm build` runs the same mapping and can add PeptideAtlas and CPTAC phosphosites
to the table.

**Build**:

```bash
# Download, map and build (UniProt sites)
hvantk reprocess uniprot-ptm:sites \
  --raw-dir data/ptm/ \
  --output data/ptm/ptm_sites.ht

# From files already in --raw-dir (a UniProt TSV and, optionally, the GTF)
hvantk reprocess uniprot-ptm:sites \
  --skip-download \
  --raw-dir data/ptm/ \
  --output data/ptm/ptm_sites.ht

# Multi-source build, with pre-downloaded files
hvantk ptm build \
  --gtf-path data/ref/Homo_sapiens.GRCh38.113.gtf.gz \
  --ptm-tsv data/ptm/uniprot-ptm-human-<YYYY-MM-DD>.tsv \
  --output-dir data/ptm/ \
  --output-ht data/ptm/ptm_sites.ht
```

With `--skip-download`, the newest `uniprot-ptm-human-<YYYY-MM-DD>.tsv` in `--raw-dir`
is mapped (`hvantk download uniprot-ptm --output-dir data/ptm/` writes one), and the GTF
is downloaded only if `Homo_sapiens.GRCh38.113.gtf.gz` is not already there. The mapping
uses GRCh38 coordinates, so `--plugin-arg reference_genome` accepts only `GRCh38`.

### Ensembl GTF (~64 MB compressed)

The same pinned-release GTF as `ensembl-gene:structure` above, used by the PTM mapper for
residue-to-genomic coordinate mapping (exon coordinates and CDS phases).
URL: https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/

**Download**: Auto-downloaded by `hvantk reprocess uniprot-ptm:sites` (into its
`--raw-dir`) and by `hvantk ptm build`. For manual download:

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
  --output cptac_phospho_brca.h5ad \
  --plugin-arg cancer_type=brca
```

**`cptac:expression`** (whole-proteome protein expression, also one AnnData per cancer
type) ships the same shape of downloader (`hvantk download cptac-expression`; also in
the built-in-downloader list above):

```bash
hvantk reprocess cptac:expression \
  --raw-dir data/cptac/ \
  --output cptac_expression_brca.h5ad \
  --plugin-arg cancer_type=brca
```

The `cancer_type` plugin-arg flows to both the download (fetch only that cancer) and the
build (select its matrix+metadata pair). `--all` on `hvantk download cptac-expression`
fetches every cancer type but does not merge them -- unlike `cptac:phospho`'s pan-cancer
TSV merge, each cancer type's expression matrix keeps its own gene set and sample
columns. A long-format TSV/CSV staged by hand still works too, via the
`{"expression","metadata"}` mapping form of `parsed_input` -- see
`hvantk/skills/cptac/expression/SKILL.md` §6-7 for the full workflow.

### PeptideAtlas phospho sites

Human phosphorylation sites from PeptideAtlas (`peptideatlas:phospho`), built from the
PeptideAtlas relational TSV dump.
URL: https://peptideatlas.org/

```bash
hvantk reprocess peptideatlas:phospho \
  --raw-dir data/peptideatlas/ \
  --output peptideatlas_phospho.parquet
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
  --output onek_genomes_variants.mt \
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

> **Caveat:** eQTLGen distributes a Z-score, not a beta/SE, and no allele frequency —
> `beta`, `se`, `af` and `maf` all come back null, so this source cannot feed
> colocalisation and only supports p-value-based overlap. Its coordinates are GRCh37
> (b37); pass `--plugin-arg reference_genome=GRCh37` to match them, or the contig
> format and the build label in the constructed variant ID disagree.

```bash
# Place cis-eQTLs_full.txt.gz under /data/eqtlgen/ then:
hvantk reprocess gtex-eqtl:eqtls \
  --raw-dir /data/eqtlgen/ \
  --output eqtl_blood.ht \
  --skip-download \
  --plugin-arg source=eqtlgen \
  --plugin-arg reference_genome=GRCh37
```

### Fang et al. (2025) pQTL data

Protein quantitative trait loci from Fang et al. (2025), covering 5 tissues (Colon, Heart, Liver, Lung, Thyroid). Space-delimited allpairs format with columns: `gene_name SNP CHR BP A1 NMISS BETA STAT P`. SE is derived as `|BETA/STAT|` (rows with `STAT = 0` are filtered out).

No direct data URL: the allpairs files are supplementary material to a preprint, not a
published dataset. Source: Fang et al. (2025), "Regulation of protein abundance in normal
human tissues", medRxiv, doi:10.1101/2025.01.10.25320181.

> **Note:** Fang pQTL data uses gene symbols. For symbol → Ensembl ID mapping, pass the HGNC lookup table built by `hvantk reprocess hgnc:lookup` via `--plugin-arg hgnc_ht=<path>`. Symbols that are not approved HGNC symbols are resolved through previous symbols, then aliases, when exactly one gene lists them; the build fails if no symbol maps.

```bash
# Build pQTL table with gene mapping
# Place Liver.allpairs_nobsGE72.txt.gz under /data/fang_pqtl/ then (the tissue name
# is inferred from the filename prefix before the first dot, so this exact naming
# matters -- see hvantk/core/utils/qtl_helpers.py):
hvantk reprocess pqtl:metrics \
  --raw-dir /data/fang_pqtl/ \
  --output pqtl_liver.ht \
  --skip-download \
  --plugin-arg source=gtex_fang \
  --plugin-arg tissue=Liver \
  --plugin-arg hgnc_ht=hgnc.ht \
  --plugin-arg p_threshold=5e-8

# Allpairs for coloc (omit p_threshold to keep all variants)
hvantk reprocess pqtl:metrics \
  --raw-dir /data/fang_pqtl/ \
  --output pqtl_allpairs_liver.ht \
  --skip-download \
  --plugin-arg source=gtex_fang \
  --plugin-arg tissue=Liver \
  --plugin-arg hgnc_ht=hgnc.ht
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
