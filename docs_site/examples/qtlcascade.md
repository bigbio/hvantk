# QTL Cascade Example

> **Heads up — build examples need refresh.** Sections that show `hvantk mktable eqtl` / `hvantk mktable pqtl` reference retired CLIs. The unified replacement is `hvantk reprocess <plugin>:<dataset>` — see the [Usage Guide](../guide/usage.md#1-build-a-dataset-with-hvantk-reprocess).

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
# Step 1: Build eQTL table from GTEx v11 significant pairs
hvantk mktable eqtl \
  --raw-input /data/gtex_v11/signif_pairs/Liver.v11.signif_pairs.parquet \
  --output-ht /data/tables/eqtl_liver.ht \
  --source gtex_v11 \
  --tissue Liver

# Step 2: Build pQTL table from Fang et al. allpairs
hvantk mktable pqtl \
  --raw-input /data/fang_pqtl/Liver_allpairs.txt.gz \
  --output-ht /data/tables/pqtl_liver.ht \
  --source gtex_fang \
  --tissue Liver \
  --gene-map-ht /data/tables/ensembl_gene.ht \
  --p-threshold 5e-8

# Step 3: Run the cascade pipeline
hvantk qtlcascade run \
  --eqtl-ht /data/tables/eqtl_liver.ht \
  --pqtl-ht /data/tables/pqtl_liver.ht \
  --constraint-ht /data/tables/gnomad_metrics.ht \
  -o /results/cascade_liver
```

## Quick Start (Python API)

```python
from hvantk.qtlcascade import (
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
# Build allpairs tables (set p-threshold to 0 to keep all variants)
hvantk mktable eqtl \
  --raw-input /data/gtex_v11/allpairs/Liver/ \
  --output-ht /data/tables/eqtl_allpairs_liver.ht \
  --source gtex_v11 --tissue Liver --p-threshold 0

hvantk mktable pqtl \
  --raw-input /data/fang_pqtl/Liver_allpairs.txt.gz \
  --output-ht /data/tables/pqtl_allpairs_liver.ht \
  --source gtex_fang --tissue Liver \
  --gene-map-ht /data/tables/ensembl_gene.ht

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

## Data Sources

| Source | Type | Format | Reference |
|--------|------|--------|-----------|
| GTEx v11 | eQTL | Parquet | GTEx Consortium |
| GTEx v8 | eQTL | TSV | GTEx Consortium |
| eQTLGen | eQTL | TSV | Vosa et al. (2021) |
| Fang et al. 2025 | pQTL | Space-delimited TSV | Fang et al. (2025) |

See [Data Sources](../guide/data-sources.md#qtl-data) for download instructions.

---

**Documentation:** [QTL Cascade Docs](../tools/qtlcascade.md)
