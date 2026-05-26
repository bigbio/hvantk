# hvantk Usage Guide

This guide covers how to build datasets and run analysis tools in hvantk.

> **Heads up:** hvantk recently retired the `mktable` and `mkmatrix` CLIs. The unified replacement is `hvantk reprocess <plugin>:<dataset>`, which runs the full Phase B pipeline (download → parse → build → drift check → save with provenance). Tutorial pages under `docs_site/guide/` and `docs_site/tools/` that still reference `mktable` / `mkmatrix` are being refreshed; see section 1 below for the current pattern.

If you haven't installed hvantk yet, see the main README for install steps.

For downloading raw data files (built-in downloaders and manual steps), see [Data Sources](data-sources.md).

## 1) Build a dataset with `hvantk reprocess`

Every in-tree provider is a plugin under `hvantk/skills/<provider>/` with a `plugin.yaml` manifest declaring its datasets, lifecycle (download / parse / build), drift probe, and artifact contract. The single CLI entry point is `hvantk reprocess`.

### Discover available datasets

```bash
# List all loaded plugins and their datasets
hvantk plugins list

# Show one plugin's manifest in detail
hvantk plugins describe clinvar
```

### Build a dataset end-to-end (download → parse → build)

```bash
# ClinVar variants
hvantk reprocess clinvar:variants \
  --raw-dir /data/clinvar/ \
  --output /out/clinvar.ht

# HGNC lookup
hvantk reprocess hgnc:lookup \
  --raw-dir /data/hgnc/ \
  --output /out/hgnc.ht
```

### Skip stages when you already have intermediates

```bash
# Already-downloaded raw file; skip the download stage
hvantk reprocess clinvar:variants \
  --skip-download \
  --raw-dir /data/clinvar/ \
  --output /out/clinvar.ht

# Pre-parsed intermediate file; skip both download and parse
hvantk reprocess gevir:metrics \
  --skip-download --skip-parse \
  --intermediate /data/gevir.tsv.bgz \
  --output /out/gevir.ht
```

### Pass plugin-specific parameters

Use `--plugin-arg KEY=VALUE` (repeatable) to forward arguments that the plugin's download / parse / build functions accept. For example, to fetch a specific CPTAC cancer type:

```bash
hvantk reprocess cptac:expression \
  --raw-dir /data/cptac/ \
  --output /out/cptac_brca.h5ad \
  --plugin-arg cancer_type=brca
```

> The deep-dive per-source tutorials under [Data Sources](data-sources.md) and the `tools/` pages still reference the retired `mktable` / `mkmatrix` CLIs; they will be rewritten in a follow-up. Use the patterns above for new work.

## 2) Ancestry Inference

Predict genetic ancestry for samples using PCA and Random Forest classification against a labeled reference panel.

### Basic Usage

```bash
# Predict ancestry using 1000 Genomes as reference
hvantk ancestry-inference \
  -q /data/my_cohort.mt \
  -r /data/1kg_phase3.mt \
  --ancestry-col super_pop \
  -o /out/ancestry_predictions.ht \
  --generate-report \
  --export-tsv
```

### With Custom Parameters

```bash
# Conservative assignment with custom filtering
hvantk ancestry-inference \
  -q /data/my_cohort.mt \
  -r /data/1kg_phase3.mt \
  --ancestry-col super_pop \
  -o /out/ancestry_predictions.ht \
  --min-af 0.05 \
  --min-call-rate 0.99 \
  --n-pcs 30 \
  --n-pcs-classify 15 \
  --min-prob 0.90 \
  --generate-report
```

### Python API

```python
import hail as hl
from hvantk.ancestry import run_ancestry_inference

# Initialize Hail
hl.init()

# Load data
query_mt = hl.read_matrix_table("my_cohort.mt")
reference_mt = hl.read_matrix_table("1kg_phase3.mt")

# Run inference
result = run_ancestry_inference(
    query_mt=query_mt,
    reference_mt=reference_mt,
    ancestry_col="super_pop",
    min_prob=0.75,
)

# Get results
predictions = result.get_predictions_df()
print(predictions['predicted_ancestry'].value_counts())

# Generate visualizations
result.generate_report("ancestry_report.html")
fig = result.plot_pca()
fig.savefig("pca_plot.png", dpi=300)

# Annotate original MT with predictions
annotated_mt = result.annotate_matrixtable(query_mt)
```

### Output Files

| File | Description |
| ---- | ----------- |
| `predictions.ht` | Hail Table with ancestry predictions |
| `predictions.tsv` | TSV export (with `--export-tsv`) |
| `ancestry_report.html` | HTML report with visualizations |
| `rf_model.pkl` | Trained model (with `--save-model`) |
| `pca_loadings.ht` | PCA loadings (with `--save-loadings`) |

**[Full Ancestry Documentation](../tools/ancestry.md)** | **[Examples](../examples/ancestry.md)**

## Expression Analysis

Inspect, summarize, and extract markers from expression MatrixTables.

```bash
# Inspect column metadata fields
hvantk expression describe -m /data/heart_sc.mt

# Collapse into gene-level summary grouped by cell type
hvantk expression summarize \
  -m /data/heart_sc.mt \
  --group-by cell_type \
  -o /out/heart_celltype_summary.ht

# Multi-field grouping with pre-filtering
hvantk expression summarize \
  -m /data/heart_sc.mt \
  --group-by cell_type --group-by region \
  --filter-by time_point=9wpc \
  --min-cells 50 \
  -o /out/heart_summary.ht

# Extract marker genes (fold-change method)
hvantk expression markers \
  -s /out/heart_celltype_summary.ht \
  --method fold_change \
  --top-n 200 \
  -o /out/heart_markers.json

# Extract markers using Wilcoxon rank-sum test
hvantk expression markers \
  -m /data/heart_sc.mt \
  --method wilcoxon \
  --group-by cell_type \
  --top-n 200 \
  -o /out/heart_wilcoxon_markers.json
```

## Prepare Custom Gene Sets

Convert plain-text gene panels into `GeneSetCollection` JSON files for use
with `hvantk psroc`, `hvantk enrichex burden`, and `hvantk enrichex overlap`.

### CLI

```bash
# Format: headerless two-column TSV (gene_set_name<TAB>gene_symbol)
hvantk genesets prepare -i panels.tsv -o panels.json

# With HGNC validation and alias resolution (recommended)
hvantk genesets prepare -i panels.tsv -o panels.json --hgnc /data/hgnc.ht

# Filter small sets and provide explicit background
hvantk genesets prepare -i panels.tsv -o panels.json \
  --hgnc /data/hgnc.ht --min-genes 5 --background bg_genes.txt

# Also export as GMT for GSEA compatibility
hvantk genesets prepare -i panels.tsv -o panels.json --export-gmt panels.gmt

# Extract gene sets from COSMIC Cancer Gene Census
hvantk genesets cosmic --ht /data/cosmic_cgc.ht -o cosmic_gene_sets.json
```

### Python API

```python
from hvantk.utils.geneset_io import parse_geneset_tsv, validate_with_hgnc
from hvantk.utils.gene_sets import load_gene_sets_from_dict

# Parse TSV
result = parse_geneset_tsv("panels.tsv")

# Optional: validate against HGNC
vr = validate_with_hgnc(result.gene_sets, "/data/hgnc.ht")

# Build and save collection
collection = load_gene_sets_from_dict(vr.gene_sets, source="prepare-geneset")
collection.save("panels.json")
```

## ClinGen Gene-Disease streamer

```python
from hvantk.data.clingen_streamer import ClinGenStreamer

streamer = ClinGenStreamer("/data/clingen/clingen_gene_disease.ht")

# High-confidence genes
definitive = streamer.get_genes_by_classification("Definitive")

# Disease keyword search
cancer_genes = streamer.get_genes_by_disease(
    ["cancer", "carcinoma", "tumor"],
    match_mode="contains",
    min_classification="Moderate",
)

# Dataset stats + summary
stats = streamer.compute_stats()
summary = streamer.classification_summary()

# Export gene sets for EnrichEx
streamer.export_for_enrichex(
    "/out/clingen_gene_sets.json",
    min_classification="Moderate",
)

# Translate output to Ensembl IDs using GeneMapper
import hail as hl
from hvantk.data.gene_mapper import GeneMapper

hgnc_ht = hl.read_table("/data/hgnc/hgnc.ht")
mapper = GeneMapper(hgnc_ht)

ensembl_ids = streamer.get_genes_by_classification(
    "Definitive",
    gene_mapper=mapper,
    output_id_type="ensembl_gene_id",
)

# Or translate to HGNC IDs
hgnc_ids = streamer.to_gene_set(
    min_classification="Moderate",
    gene_mapper=mapper,
    output_id_type="hgnc_id",
)
```

## File Format Conversion

Hail supports standard gzip (`.gz`) and uncompressed files but processes them single-threaded. Block gzip (BGZF) `.bgz` files enable parallel import and are strongly recommended for large datasets. Convert with `hvantk utils convert-bgz`:

```bash
# Default: replaces .gz extension with .bgz
hvantk utils convert-bgz input.tsv.gz

# Custom output path and thread count
hvantk utils convert-bgz input.tsv.gz -o output.tsv.bgz --threads 4
```

The command auto-detects whether the file is already BGZF and skips conversion if so.

## Tips & troubleshooting

- Use `--overwrite` to replace an existing output. Without it, builders abort if the output exists.
- For JSON vs YAML: JSON works out of the box; YAML recipes require `PyYAML` installed.
- For UCSC, gene labels may be pipe-delimited (e.g., A|B); `--split-gene-field` defaults to true.
- MatrixTables typically store sample/cell metadata under `mt.col_key` and cols metadata; inspect with `mt.describe()` in Python or logs from CLI.
- **gzip vs BGZF**: Hail reads standard gzip files single-threaded, which is significantly slower for large files. Convert to BGZF with `--auto-convert-bgz` or `hvantk utils convert-bgz` for parallel import.

## See also

- [Architecture](../architecture.md) – system design and extension points
- [Contributing](../contributing.md) – development workflow and contribution guidelines
- [Recipe Examples](https://github.com/bigbio/hvantk/tree/main/examples/recipes/) – ready-to-edit recipe templates
