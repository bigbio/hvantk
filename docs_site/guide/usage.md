# hvantk Usage Guide

This guide covers how to build datasets and run analysis tools in hvantk.

> **Migration note:** hvantk has retired the `mktable` and `mkmatrix` CLIs. The unified replacement is `hvantk reprocess <plugin>:<dataset>`, which runs the full Phase B pipeline (download → parse → build → drift check → save with provenance). See section 1 below for the current pattern.

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

CPTAC phosphoproteomics builds one AnnData per cancer type — the single
`cancer_type` flows to both the download and the build:

```bash
hvantk reprocess cptac:phospho \
  --raw-dir /data/cptac/ \
  --output /out/cptac_phospho_brca.h5ad \
  --plugin-arg cancer_type=brca
```

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
from hvantk.algorithms.ancestry import run_ancestry_inference

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

Inspect, summarize, and extract markers from expression AnnData (`.h5ad`) files.

```bash
# Inspect observation metadata fields
hvantk expression describe -m /data/heart_sc.h5ad

# Collapse into a per-group × per-gene summary AnnData grouped by cell type
hvantk expression summarize \
  -m /data/heart_sc.h5ad \
  --group-by cell_type \
  -o /out/heart_celltype_summary.h5ad

# Multi-field grouping with pre-filtering
hvantk expression summarize \
  -m /data/heart_sc.h5ad \
  --group-by cell_type --group-by region \
  --filter-by time_point=9wpc \
  --min-cells 50 \
  -o /out/heart_summary.h5ad

# Extract marker genes (Wilcoxon rank-sum test)
hvantk expression markers \
  -m /data/heart_sc.h5ad \
  --method wilcoxon \
  --group-by cell_type \
  --top-n 200 \
  -o /out/heart_markers.json

# Extract markers using a t-test, pre-filtered to one region
hvantk expression markers \
  -m /data/heart_sc.h5ad \
  --method t-test \
  --group-by cell_type \
  --filter-by region=LV \
  --top-n 200 \
  -o /out/heart_LV_markers.json
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
from hvantk.core.utils.geneset_io import parse_geneset_tsv, validate_with_catalog
from hvantk.core.utils.gene_sets import load_gene_sets_from_dict
from hvantk.skills.hgnc.streamers import HGNCGeneCatalogStreamer

# Parse TSV
result = parse_geneset_tsv("panels.tsv")

# Optional: validate against the HGNC gene catalog (resolves aliases)
catalog = HGNCGeneCatalogStreamer.from_path("/data/hgnc.ht")
vr = validate_with_catalog(result.gene_sets, catalog)

# Build and save collection
collection = load_gene_sets_from_dict(vr.gene_sets, source="prepare-geneset")
collection.save("panels.json")
```

## ClinGen Gene-Disease streamer

```python
from hvantk.skills.clingen.streamers import ClinGenGeneDiseaseTableStreamer

streamer = ClinGenGeneDiseaseTableStreamer("/data/clingen/clingen_gene_disease.ht")

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

# Translate output to Ensembl IDs using the HGNC gene catalog
from hvantk.skills.hgnc.streamers import HGNCGeneCatalogStreamer

catalog = HGNCGeneCatalogStreamer.from_path("/data/hgnc/hgnc.ht")

ensembl_ids = streamer.get_genes_by_classification(
    min_classification="Definitive",
    gene_catalog=catalog,
    output_id_type="ensembl_gene_id",
)

# Or translate to HGNC IDs
hgnc_ids = streamer.to_gene_set(
    min_classification="Moderate",
    gene_catalog=catalog,
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

## Gene-level annotation matrices

`hvantk annotate` turns built datasets into one gene × feature matrix. It is three
commands run in order, because each stage is independently re-runnable and the
expensive middle stage is per-source:

| Stage | Command | Produces |
| --- | --- | --- |
| 1. Spine | `annotate spine` | the gene table every annotation joins onto |
| 2. Prepare | `annotate prepare` (once per axis) | one source mapped onto the spine's `gene_id` |
| 3. Compose | `annotate compose` | the joined matrix + a manifest JSON |

### 1. Build the spine

The spine fixes which genes exist and what `gene_id` means for the whole matrix.

```bash
hvantk annotate spine \
  --gene-table ensembl_structure.ht \
  --hgnc hgnc_lookup.ht \
  --output spine.ht
```

`--biotype` defaults to `protein_coding`; pass `--biotype all` to keep every biotype.

### 2. Declare the sources

A feature spec names each axis, the dataset it comes from, the key to join on, and the
columns to keep:

```yaml
name: chd-features
layer1:
  - {axis: constraint, source: gnomad-metrics:metrics, key: gene_id, columns: [mis_z]}
  - {axis: gevir,      source: gevir:table,            key: gene_id, columns: [gevir_pct]}
```

Each `axis` label must be unique — it is the identity the later stages join on. When an
entry's `key` is `hgnc_id` or `symbol` rather than `gene_id`, `prepare` needs
`--hgnc` so it can resolve the mapping.

### 3. Prepare each axis, then compose

`prepare` runs once per axis and is the stage worth parallelising — each invocation is
independent, so they can be separate cluster jobs:

```bash
hvantk annotate prepare --spec features.yaml --axis constraint \
  --input gnomad_metrics.ht --spine spine.ht --output constraint.ht

hvantk annotate prepare --spec features.yaml --axis gevir \
  --input gevir.ht --spine spine.ht --output gevir.ht
```

`compose` left-joins every prepared axis onto the spine. It requires one `--prepared`
per axis in the spec and fails if any is missing, so a partial matrix cannot be produced
silently:

```bash
hvantk annotate compose --spec features.yaml --spine spine.ht \
  --prepared constraint=constraint.ht \
  --prepared gevir=gevir.ht \
  --output features.ht
```

The manifest JSON lands next to the output (`features.ht.manifest.json` unless
`--manifest` says otherwise) and records what went into the matrix.

Scheduling is deliberately external: the commands do their work in-process and hvantk
never imports a scheduler, so an sbatch wrapper submits them.

## External cohorts

`hvantk cohort` brings a cohort's own gene-level results alongside the annotation
matrix. A cohort is declared by a manifest, not by flags:

```yaml
name: my-cohort
key: symbol              # or gene_id / hgnc_id
key_column: gene         # the column in `table` holding that key
table: cohort.tsv
prior:
  column: minp
  direction: lower_is_better
cohort_axes:
  - axis: architecture
    columns: [n_case_var, conc, driver_af]
```

```bash
# Check the manifest against its table (and labels) before anything reads it
hvantk cohort validate --manifest cohort.yaml

# Originate a gene-level prior with a Fisher-exact burden test
hvantk cohort burden --help

# Join the cohort's declared columns onto the Layer-1 matrix
hvantk cohort attach --help
```

`validate` first is the point: it fails on a key column that is not present, or a
declared axis column missing from the table, rather than letting a silently-empty join
propagate into the matrix.

The same manifest is what `hvantk rerank` consumes as its starting ranking — see
[examples/rerank](https://github.com/bigbio/hvantk/tree/main/examples/rerank) for a
complete runnable configuration.

## Tips & troubleshooting

- Use `--overwrite` to replace an existing output. Without it, builders abort if the output exists.
- For UCSC, gene labels may be pipe-delimited (e.g., A|B); `--split-gene-field` defaults to true.
- Expression AnnData stores per-cell/sample metadata in `adata.obs`; inspect it with `hvantk expression describe -m <file.h5ad>` or `adata.obs` in Python.
- **gzip vs BGZF**: Hail reads standard gzip files single-threaded, which is significantly slower for large files. Pre-convert with `hvantk utils convert-bgz <file.gz>` for parallel import before running `hvantk reprocess`.

## See also

- [Architecture](../architecture.md) – system design and extension points
- [Contributing](../contributing.md) – development workflow and contribution guidelines
- [Data Sources](data-sources.md) – available datasets and how to acquire them
