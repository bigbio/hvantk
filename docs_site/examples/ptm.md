# PTM (Post-Translational Modification) Variant Classification Example

This page walks through the end-to-end PTM pipeline: building a PTM sites table, annotating variants, running landscape and population analyses, and generating a report.

## Overview

The PTM pipeline maps UniProt post-translational modification sites to GRCh38 genomic coordinates, cross-references them with ClinVar and gnomAD variants, and quantifies the relationship between PTM sites and variant pathogenicity.

## Prerequisites

```bash
# Activate environment
eval "$(poetry env activate)"
```

### Required Data

| Input | Source | How to obtain |
|-------|--------|---------------|
| ClinVar Hail Table | NCBI | `hvantk download clinvar` then `hvantk mktable clinvar` |
| gnomAD Hail Table | gnomAD | Manual download (see [Data Sources](../guide/data-sources.md)) |
| Ensembl GTF | Ensembl | Auto-downloaded by `hvantk ptm build` |
| UniProt PTM TSV | UniProt | Auto-downloaded by `hvantk ptm build` |

## Step 1: Build PTM Sites Table

Download PTM data, map to genomic coordinates, and build a Hail Table.

```bash
# Automatic download and build
hvantk ptm build \
  --output-dir data/ptm/ \
  --output-ht data/ptm/ptm_sites.ht

# Or with pre-downloaded inputs
hvantk ptm build \
  --gtf-path data/ref/Homo_sapiens.GRCh38.113.gtf.gz \
  --ptm-tsv data/ptm/uniprot-ptm-human.tsv \
  --output-dir data/ptm/ \
  --output-ht data/ptm/ptm_sites.ht
```

**Expected output:**
```
Mapping complete: 55422/56120 mapped (98.7%), 698 failed
Resolution: {'ensembl_canonical': 48231, 'longest_cds': 7191}
Mapped TSV: data/ptm/ptm_mapped.tsv
Hail Table: data/ptm/ptm_sites.ht
```

## Step 2: Annotate Variants

Cross-reference ClinVar variants with PTM sites.

```bash
hvantk ptm annotate \
  --variants-ht data/clinvar.ht \
  --ptm-ht data/ptm/ptm_sites.ht \
  -o data/clinvar_ptm.ht
```

**Expected output:**
```
Annotated 200,000 variants:
  PTM site:  1,234
  Proximal:  5,678
Output: data/clinvar_ptm.ht
```

## Step 3: Landscape Analysis (Q1)

PTM-variant overlap and enrichment analysis.

```bash
hvantk ptm landscape \
  --clinvar-ht data/clinvar.ht \
  --ptm-ht data/ptm/ptm_sites.ht \
  -o results/landscape/ \
  --save-plots
```

**Expected output:**
```
PTM-Variant Landscape:
  Variants: 200,000 total, 50,000 P/LP, 30,000 B/LB
  P/LP at PTM site: 1,234
  P/LP proximal: 5,678
  Enrichment: OR=2.15, p=1.23e-45
  By PTM category (P/LP):
    Phosphorylation: 456
    Ubiquitination: 234
    Acetylation: 123
Plots saved to results/landscape/
```

**Output files:**
- `results/landscape/landscape_summary.json`
- `results/landscape/landscape_summary.png`
- `results/landscape/overlap_by_category.png`
- `results/landscape/distance_distribution.png`

## Step 4: Predictor Evaluation (Q2 - Composed Workflow)

Export PTM-stratified variant lists, then run PSROC independently on each stratum.

```bash
# Export strata
hvantk ptm export-strata --annotated-ht data/clinvar_ptm.ht -o strata/

# Run PSROC on PTM variants
hvantk psroc \
  --variants strata/ptm_variants.txt \
  --clinvar-ht data/clinvar.ht \
  --dbnsfp-ht data/dbnsfp.ht \
  --scores "CADD_phred,REVEL_score,MetaLR_score" \
  --output-dir results/psroc_ptm/

# Run PSROC on non-PTM variants
hvantk psroc \
  --variants strata/non_ptm_variants.txt \
  --clinvar-ht data/clinvar.ht \
  --dbnsfp-ht data/dbnsfp.ht \
  --scores "CADD_phred,REVEL_score,MetaLR_score" \
  --output-dir results/psroc_non_ptm/
```

Compare AUC values between strata to determine whether predictor performance differs at PTM sites.

## Step 5: Population Analysis (Q3)

Compare allele frequency distributions at PTM sites vs non-PTM coding positions in gnomAD.

```bash
hvantk ptm population \
  --gnomad-ht data/gnomad.ht \
  --ptm-ht data/ptm/ptm_sites.ht \
  -o results/population/ \
  --save-plots

# Optionally include CCR scores
hvantk ptm population \
  --gnomad-ht data/gnomad.ht \
  --ptm-ht data/ptm/ptm_sites.ht \
  --ccr-ht data/ccr.ht \
  -o results/population/ \
  --save-plots
```

**Expected output:**
```
PTM Population Analysis:
  Total variants: 5,000,000
  At PTM site: 12,345 (mean AF=1.23e-04)
  Proximal: 45,678 (mean AF=2.34e-04)
  Non-PTM: 4,941,977 (mean AF=5.67e-04)
  PTM sites with zero AF: 8,901
```

## Step 6: Generate Report

Combine landscape and population results into a single HTML report.

```bash
hvantk ptm report -o results/ptm_report.html \
  --landscape-json results/landscape/landscape_summary.json \
  --population-json results/population/population_summary.json \
  --title "PTM-Variant Classification Report" \
  --description "ClinVar and gnomAD analysis of UniProt PTM sites"
```

The report includes summary cards, embedded plots, per-category tables, and a methods section.

## Python API

All steps can also be run programmatically:

```python
import hail as hl
from hvantk.ptm import (
    PTMBuildConfig,
    ptm_build_pipeline,
    annotate_variants_with_ptm,
    ptm_landscape,
    ptm_population,
    export_ptm_strata,
)
from hvantk.ptm.report import generate_report

# Build
config = PTMBuildConfig(output_dir="data/ptm/", output_ht="data/ptm/ptm_sites.ht")
build_result = ptm_build_pipeline(config)

# Annotate
clinvar = hl.read_table("data/clinvar.ht")
ptm = hl.read_table("data/ptm/ptm_sites.ht")
annotated = annotate_variants_with_ptm(clinvar, ptm)
annotated = annotated.checkpoint("data/clinvar_ptm.ht")

# Landscape (Q1)
landscape = ptm_landscape(clinvar, ptm, "results/landscape/")

# Export strata (Q2)
strata = export_ptm_strata(annotated, "strata/")

# Population (Q3)
gnomad = hl.read_table("data/gnomad.ht")
population = ptm_population(gnomad, ptm, "results/population/")

# Report
generate_report(
    "results/ptm_report.html",
    landscape_result=landscape,
    population_result=population,
)
```

## Documentation

- [PTM Documentation](../tools/ptm.md)
- [PSROC Documentation](../tools/psroc.md) (for Q2 predictor evaluation)
- [Data Sources](../guide/data-sources.md)
- [Usage Guide](../guide/usage.md)
