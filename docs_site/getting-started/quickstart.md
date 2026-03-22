# Quick Start

This page provides an overview of the core workflows available in hvantk. Each section includes a minimal CLI example to get you started — follow the links for full documentation.

## HGC: Joint Genotyping Pipeline

High-performance joint genotyping for large cohorts. Combines thousands of GVCF files with integrated QC.

```bash
# End-to-end pipeline
hvantk hgc pipeline -i /data/gvcfs -o /output

# Or run individual steps
hvantk hgc gvcf-combine -g /data/gvcfs -o cohort.vds
hvantk hgc vds2mt -i cohort.vds -o cohort.mt
hvantk hgc qc-report -i cohort.mt -o qc_report.html
```

[Full HGC Documentation](../tools/hgc.md){ .md-button }

## PSROC: Variant Score Evaluation

Evaluate pathogenicity prediction scores (CADD, REVEL, MetaLR) using ClinVar truth labels. Generate ROC curves and performance metrics.

```bash
hvantk psroc \
  --genes-file genes.txt \
  --clinvar-ht clinvar.ht \
  --dbnsfp-ht dbnsfp.ht \
  --scores "CADD_phred,REVEL_score" \
  --output-dir results/
```

[PSROC Documentation](../tools/psroc.md){ .md-button } [Example](../examples/psroc.md){ .md-button }

## EnrichEx: Gene Set Enrichment

Test gene set enrichment using overlap analysis (Fisher's exact test) and case-control burden testing (rare variant regression).

```bash
# Overlap enrichment
hvantk enrichex overlap \
  -g gwas_genes.txt \
  -s gene_sets.json \
  -o overlap_results.tsv \
  --generate-report

# Burden testing
hvantk enrichex burden \
  -m cohort.mt \
  -p phenotypes.ht \
  -s gene_sets.json \
  -o burden_results.tsv \
  --generate-report
```

[EnrichEx Documentation](../tools/enrichex.md){ .md-button } [Example](../examples/enrichex.md){ .md-button }

## Ancestry Inference

Predict genetic ancestry for samples using PCA and Random Forest classification against a labeled reference panel (e.g., 1000 Genomes).

```bash
hvantk ancestry-inference \
  -q cohort.mt \
  -r 1kg_reference.mt \
  --ancestry-col super_pop \
  -o ancestry_predictions.ht \
  --generate-report
```

[Ancestry Documentation](../tools/ancestry.md){ .md-button } [Example](../examples/ancestry.md){ .md-button }

## Annotation Tables

Create Hail Tables from public databases (ClinVar, gnomAD, Ensembl).

```bash
# Single table
hvantk mktable clinvar --raw-input clinvar.vcf.bgz --output-ht clinvar.ht

# Batch processing
hvantk mktable-batch --recipe tables_recipe.json
```

[Tables Guide](../guide/usage.md){ .md-button }

## Expression Matrices

Build Hail MatrixTables from bulk and single-cell expression data.

```bash
# UCSC Cell Browser data
hvantk mkmatrix ucsc -e expr.tsv.bgz -m metadata.tsv -o ucsc.mt

# Batch processing
hvantk mkmatrix-batch --recipe matrices_recipe.json
```

[Expression Guide](../guide/usage.md){ .md-button }

## Expression Analysis

Inspect, summarize, and extract marker genes from expression MatrixTables.

```bash
# Inspect metadata
hvantk expression describe -m ucsc.mt

# Summarize by cell type
hvantk expression summarize -m ucsc.mt --group-by cell_type -o summary.ht

# Extract marker genes
hvantk expression markers -s summary.ht --method fold_change -o markers.json
```

[Expression Guide](../guide/usage.md){ .md-button }

## File Format Conversion

Convert standard gzip files to BGZF for Hail parallel import:

```bash
hvantk utils convert-bgz input.tsv.gz -o output.tsv.bgz --threads 4
```

Or use `--auto-convert-bgz` in supported commands to convert on-the-fly.

## Data Downloaders

Download curated datasets from public repositories.

```bash
hvantk download ucsc --dataset adultPancreas --output-dir data/ucsc
```

[Data Sources](../guide/data-sources.md){ .md-button }

## Full Quick Start Example

```bash
# Download and process expression data
hvantk download ucsc --dataset adultPancreas --output-dir data/ucsc
hvantk mkmatrix ucsc -e data/ucsc/exprMatrix.tsv.bgz -m data/ucsc/meta.tsv -o data/ucsc/adultPancreas.mt

# If expression matrix is plain gzip (.gz), use --auto-convert-bgz
hvantk mkmatrix ucsc -e data/ucsc/exprMatrix.tsv.gz -m data/ucsc/meta.tsv -o data/ucsc/adultPancreas.mt --auto-convert-bgz

# Build annotation tables
hvantk mktable clinvar --raw-input clinvar.vcf.bgz --output-ht clinvar.ht --ref-genome GRCh38

# Or use batch processing with recipes
hvantk mktable-batch --recipe recipe.json
```
