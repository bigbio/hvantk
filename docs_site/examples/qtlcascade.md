# QTL Cascade Example

This page demonstrates the QTL cascade pipeline for tracing variant effects from DNA to RNA (eQTL) to protein (pQTL).

## Overview

The QTL cascade pipeline:
1. Builds eQTL and pQTL Hail Tables from summary statistics
2. Outer-joins them on `(locus, alleles, gene_id)` to classify variant–gene pairs
3. Optionally runs colocalization ABF to distinguish true cascades from LD artifacts
4. Generates gene-level summaries with constraint and disease overlays
5. Produces plots and an HTML report

## Quick Start (CLI)

```bash
# Step 1: Build eQTL table from GTEx v11 significant pairs.
# Place Liver.v11.signif_pairs.parquet under /data/gtex_v11/signif_pairs/.
hvantk reprocess gtex-eqtl:eqtls \
  --raw-dir /data/gtex_v11/signif_pairs/ \
  --output /data/tables/eqtl_liver.ht \
  --skip-download \
  --plugin-arg source=gtex_v11 \
  --plugin-arg tissue=Liver

# Step 2: Build pQTL table from Fang et al. allpairs.
# Place Liver_allpairs.txt.gz under /data/fang_pqtl/.
hvantk reprocess pqtl:metrics \
  --raw-dir /data/fang_pqtl/ \
  --output /data/tables/pqtl_liver.ht \
  --skip-download \
  --plugin-arg source=gtex_fang \
  --plugin-arg tissue=Liver \
  --plugin-arg hgnc_ht=/data/tables/ensembl_gene.ht \
  --plugin-arg p_threshold=5e-8

# Step 3: Run the cascade pipeline
hvantk qtlcascade run \
  --eqtl-ht /data/tables/eqtl_liver.ht \
  --pqtl-ht /data/tables/pqtl_liver.ht \
  --constraint-ht /data/tables/gnomad_metrics.ht \
  -o /results/cascade_liver
```

## Quick Start (Python API)

```python
from hvantk.algorithms.qtlcascade import (
    build_cascade,
    build_cascade_gene_summary,
    CascadeConfig,
    CascadePipeline,
)

# Option A: Use individual functions
ht = build_cascade(
    eqtl_ht_path="/data/tables/eqtl_liver.ht",
    pqtl_ht_path="/data/tables/pqtl_liver.ht",
    output_path="/results/cascade.ht",
    tissue="Liver",
)
print(f"Cascade pairs: {ht.count()}")

gene_ht = build_cascade_gene_summary(
    cascade_ht_path="/results/cascade.ht",
    output_path="/results/gene_summary.ht",
    constraint_ht_path="/data/tables/gnomad_metrics.ht",
)
print(f"Cascade genes: {gene_ht.count()}")

# Option B: Use the pipeline
config = CascadeConfig(
    eqtl_ht="/data/tables/eqtl_liver.ht",
    pqtl_ht="/data/tables/pqtl_liver.ht",
    constraint_ht="/data/tables/gnomad_metrics.ht",
    output_dir="/results/cascade_liver",
)
pipeline = CascadePipeline(config)
result = pipeline.run()
print(f"Class counts: {result.class_counts}")
```

## With Colocalization

To distinguish true signal propagation from LD artifacts, provide allpairs tables for coloc:

```bash
# Build allpairs tables (set p_threshold to 0 to keep all variants)
hvantk reprocess gtex-eqtl:eqtls \
  --raw-dir /data/gtex_v11/allpairs/Liver/ \
  --output /data/tables/eqtl_allpairs_liver.ht \
  --skip-download \
  --plugin-arg source=gtex_v11 \
  --plugin-arg tissue=Liver \
  --plugin-arg p_threshold=0

hvantk reprocess pqtl:metrics \
  --raw-dir /data/fang_pqtl/ \
  --output /data/tables/pqtl_allpairs_liver.ht \
  --skip-download \
  --plugin-arg source=gtex_fang \
  --plugin-arg tissue=Liver \
  --plugin-arg hgnc_ht=/data/tables/ensembl_gene.ht

# Run pipeline with coloc
hvantk qtlcascade run \
  --eqtl-ht /data/tables/eqtl_liver.ht \
  --pqtl-ht /data/tables/pqtl_liver.ht \
  --eqtl-allpairs /data/tables/eqtl_allpairs_liver.ht \
  --pqtl-allpairs /data/tables/pqtl_allpairs_liver.ht \
  --constraint-ht /data/tables/gnomad_metrics.ht \
  --disease-genes-ht /data/tables/clingen.ht \
  -o /results/cascade_liver_coloc
```

## Multi-Tissue Analysis

Run across Fang et al. (2025) tissues with cross-tissue comparison:

```bash
hvantk qtlcascade run \
  --eqtl-ht /data/tables/eqtl.ht \
  --pqtl-ht /data/tables/pqtl.ht \
  --eqtl-allpairs /data/tables/eqtl_allpairs.ht \
  --pqtl-allpairs /data/tables/pqtl_allpairs.ht \
  --tissues "Liver,Heart,Lung,Colon,Thyroid" \
  --constraint-ht /data/tables/gnomad_metrics.ht \
  -o /results/cascade_multi
```

```python
# Python API for multi-tissue
config = CascadeConfig(
    eqtl_ht="/data/tables/eqtl.ht",
    pqtl_ht="/data/tables/pqtl.ht",
    eqtl_allpairs_ht="/data/tables/eqtl_allpairs.ht",
    pqtl_allpairs_ht="/data/tables/pqtl_allpairs.ht",
    tissues=["Liver", "Heart", "Lung", "Colon", "Thyroid"],
    constraint_ht="/data/tables/gnomad_metrics.ht",
    output_dir="/results/cascade_multi",
)
pipeline = CascadePipeline(config)
results = pipeline.run_collection()

for tissue, res in results.items():
    n_coloc = 0
    if res.coloc_df is not None:
        n_coloc = (res.coloc_df["H4"] > 0.8).sum()
    print(f"{tissue}: {res.n_cascade_genes} genes, {n_coloc} colocalised")
```

## Expected Outputs

```
results/cascade_liver/
├── per_tissue/
│   └── liver/
│       ├── cascade.ht/          # Variant-gene cascade pairs
│       ├── gene_summary.ht/     # Gene-level aggregation
│       ├── gene_summary.tsv     # TSV export
│       └── coloc_results.tsv    # H0-H4 posteriors (if coloc run)
├── plots/
│   ├── liver_cascade_classes.png
│   └── liver_coloc_posteriors.png
└── qtlcascade_report.html
```

## Standalone Coloc

Run colocalization independently on a set of genes:

```bash
# Extract cascade genes with both eQTL + pQTL evidence
# (e.g., from gene_summary.tsv where has_complete_cascade = true)
awk -F'\t' 'NR>1 && $7=="true" {print $1}' gene_summary.tsv > cascade_genes.txt

# Run coloc
hvantk qtlcascade coloc \
  --eqtl-allpairs /data/tables/eqtl_allpairs_liver.ht \
  --pqtl-allpairs /data/tables/pqtl_allpairs_liver.ht \
  --cascade-genes cascade_genes.txt \
  --tissue Liver \
  --window-kb 500 \
  -o coloc_results.tsv
```

## GWAS → Effector Colocalization (`gwas-coloc`)

`hvantk qtlcascade gwas-coloc` colocalizes a **GWAS** locus against cis-**eQTL** to rank the likely
effector gene(s), then (by default) confirms the lead effector with **SuSiE-RSS + `coloc.susie`**
fine-mapping. All summary statistics stream via **remote tabix** — no bulk downloads.

The fine-mapping step is what separates a genuine colocalization (`CONFIRMED`) from a
single-variant-ABF artifact (`REFUTED`). The example below shows both, anchored on a
literature-supported positive control.

### Positive control — atrial fibrillation → MYOZ1 (CONFIRMED)

`MYOZ1` (10q22) is a known cardiac effector for atrial fibrillation. Confirm it end-to-end
(FinnGen `I9_AF` × GTEx heart atrial-appendage eQTL `QTD000251`):

```bash
hvantk qtlcascade gwas-coloc \
  --endpoint I9_AF --chrom 10 --lead 73600000 \
  --eqtl QTD000251 \
  --gene ENSG00000177791 \
  --gwas-n 261395 --eqtl-n 372 \
  -o results/myoz1
```

Expected output (values are illustrative; remote data may shift slightly between releases):

```text
Region chr10:73100000-74100000  |  GWAS min-p 2.7e-14  |  51 genes tested
Top effector: ENSG00000177791 (PP4=0.812)
Fine-map: credible sets GWAS=1 eQTL=1; coloc.susie PP4=0.649
VERDICT: CONFIRMED (ABF + fine-mapping agree)
Report: results/myoz1/report_I9_AF_10_73600000.json
```

> `ENSG00000177791` is `MYOZ1`. The tool reports Ensembl gene IDs; HGNC symbol mapping is a planned
> enhancement.

### Fast tier — ABF only (offline-friendly)

The ABF ranking is pure-Python (only the core deps). Skip fine-mapping with `--no-fine-map`
(then `--gwas-n`/`--eqtl-n` are not required):

```bash
hvantk qtlcascade gwas-coloc \
  --endpoint I9_AF --chrom 10 --lead 73600000 \
  --eqtl QTD000251 --no-fine-map \
  -o results/myoz1_abf
# Top effector: ENSG00000177791 (PP4=0.812)
# VERDICT: SUGGESTIVE (ABF only; fine-mapping not run)
```

### Contrast — a CHD locus that does NOT survive fine-mapping (REFUTED)

The septal-defect 17q21 locus has an *identical-looking* ABF signal (PP4 ≈ 0.81 for `NSF`) that
**fails** fine-mapping — neither trait yields a credible set, so the ABF hit is a single-variant
artifact:

```bash
hvantk qtlcascade gwas-coloc \
  --endpoint Q17_SEPTA_DEFEC --chrom 17 --lead 46890164 \
  --eqtl QTD000136 \
  --gene ENSG00000073969 \
  --gwas-n 412181 --eqtl-n 213 \
  -o results/nsf_17q21
# Top effector: ENSG00000073969 (PP4=0.81)
# Fine-map: credible sets GWAS=0 eQTL=0; coloc.susie PP4=0.0
# VERDICT: REFUTED (no fine-mappable signal — single-variant-ABF artifact)
```

This `CONFIRMED`-vs-`REFUTED` contrast is why the fine-mapping layer matters: single-variant ABF
over-calls when a strong GWAS meets a weak eQTL.

### Prerequisites

| Tier | Requirements |
|------|--------------|
| ABF only (`--no-fine-map`) | core hvantk deps (pysam, numpy, pandas) + network |
| Fine-mapping (default) | **same core deps** — fine-mapping is pure-Python SuSiE-RSS + `coloc.susie` (no `R` / `bcftools` / `curl`). Streams a 1000G GRCh38 LD reference for the chosen `--superpop` (default `EUR`), cached under `--ld-cache-dir` (default: `$HVANTK_LD_CACHE`, else `~/.cache/hvantk/1kg`); or pass a local `--ld-vcf` for an **offline LD reference** (the GWAS/eQTL summary statistics are still streamed) |

Data sources (all remote, no downloads): FinnGen R10 GWAS, eQTL Catalogue (GTEx) cis-eQTL, and
1000 Genomes high-coverage GRCh38 for the LD reference.

### Python API

```python
from hvantk.algorithms.qtlcascade.gwas_pipeline import (
    GwasColocConfig, run_gwas_coloc_pipeline,
)

config = GwasColocConfig(
    endpoint="I9_AF", chrom="10", lead=73600000,
    eqtl_dataset="QTD000251", gene_of_interest="ENSG00000177791",
    gwas_N=261395, eqtl_N=372, fine_map=True,
    output_dir="results/myoz1",
)
report = run_gwas_coloc_pipeline(config)
print(report["verdict"])             # CONFIRMED (ABF + fine-mapping agree)
print(report["results"]["top_PP4"])  # ~0.81
```

## Data Sources

| Source | Type | Format | Reference |
|--------|------|--------|-----------|
| GTEx v11 | eQTL | Parquet | GTEx Consortium |
| GTEx v8 | eQTL | TSV | GTEx Consortium |
| eQTLGen | eQTL | TSV | Vosa et al. (2021) |
| Fang et al. 2025 | pQTL | Space-delimited TSV | Fang et al. (2025) |
| FinnGen R10 | GWAS | Remote-tabix (bgzip+tbi) | FinnGen (2023) |
| eQTL Catalogue (GTEx) | cis-eQTL | Remote-tabix (bgzip+tbi) | Kerimov et al. (2021) |
| 1000 Genomes (GRCh38) | LD reference | VCF (remote-tabix) | 1000G / NYGC (2020) |

See [Data Sources](../guide/data-sources.md#qtl-data) for download instructions.

---

**Documentation:** [QTL Cascade Docs](../tools/qtlcascade.md)
