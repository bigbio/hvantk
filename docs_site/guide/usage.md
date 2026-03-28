# hvantk Usage Guide

This guide shows practical, copy-pasteable examples to build Hail Tables (HT) and MatrixTables (MT) from explicit raw files and from recipes (JSON/YAML).

If you haven't installed hvantk yet, see the main README for install steps.

For downloading raw data files (built-in downloaders and manual steps), see [Data Sources](data-sources.md).

## 1) Build a single annotation Table (HT)

Build one table at a time with explicit inputs and options.

- ClinVar (VCF → HT keyed by [locus, alleles])

```bash
hvantk mktable clinvar \
  --raw-input /data/clinvar_2024.vcf.bgz \
  --output-ht /out/clinvar.ht \
  --ref-genome GRCh38 \
  --overwrite
```

- Interactome (BED intervals → HT keyed by interval)

```bash
hvantk mktable interactome \
  --raw-input /data/insider.bed.bgz \
  --output-ht /out/interactome.ht
```

- GeVIR (TSV keyed by gene_id)

```bash
hvantk mktable gevir \
  --raw-input /data/gevir.tsv.bgz \
  --output-ht /out/gevir.ht \
  --fields oe_syn_upper,oe_mis_upper
```

- gnomAD constraint metrics (TSV keyed by gene_id)

```bash
hvantk mktable gnomad-metrics \
  --raw-input /data/gnomad.tsv.bgz \
  --output-ht /out/gnomad.ht
```

- dbNSFP variant annotations (TSV keyed by locus, alleles)

```bash
hvantk mktable dbnsfp \
  --raw-input /data/dbNSFP4_variant.bgz \
  --output-ht /out/dbnsfp.ht \
  --ref-genome GRCh38 \
  --auto-convert-bgz
```

> **Tip:** If your dbNSFP file is standard gzip (`.gz`) rather than BGZF, add `--auto-convert-bgz` to automatically convert it before import.

- Ensembl gene annotations (Biomart TSV keyed by gene_id)

```bash
hvantk mktable ensembl-gene \
  --raw-input /data/biomart.tsv.bgz \
  --output-ht /out/ensembl.ht \
  --no-canonical
```

- ClinGen Gene-Disease validity (CSV keyed by gene_id + disease_id)

```bash
hvantk mktable clingen-gene-disease \
  --raw-input /data/clingen/Clingen-Gene-Disease-Summary-2026-03-22.csv \
  --output-ht /out/clingen_gene_disease.ht
```

- GenCC submissions (CSV keyed by gene_id + disease_id)

```bash
hvantk mktable gencc-submissions \
  --raw-input /data/gencc/gencc-submissions.csv \
  --output-ht /out/gencc_submissions.ht
```

- COSMIC Cancer Gene Census (TSV keyed by gene_id)

```bash
hvantk mktable cosmic-cgc \
  --raw-input /data/cosmic/cancer_gene_census.tsv \
  --output-ht /out/cosmic_cgc.ht \
  --mutation-context somatic
```

- HGNC gene nomenclature (TSV keyed by hgnc_id)

```bash
hvantk mktable hgnc \
  --raw-input /data/hgnc/hgnc_complete_set.tsv \
  --output-ht /out/hgnc.ht
```

- PTM sites (UniProt PTM → genomic coordinates, keyed by locus)

```bash
hvantk ptm build \
  --output-dir /data/ptm/ \
  --output-ht /out/ptm_sites.ht
```

> **Note:** The PTM build command downloads Ensembl GTF and UniProt PTM data automatically. Use `--gtf-path` and `--ptm-tsv` to provide pre-downloaded files.

## 2) Batch-create Tables (HT) from a recipe

Use a recipe to build many tables at once. JSON and YAML are both supported (YAML requires PyYAML installed).

Example JSON recipe (save as examples/recipes/tables.example.json):

```json
{
  "tables": [
    {
      "name": "clinvar",
      "input": "/data/clinvar_2024.vcf.bgz",
      "output": "/out/clinvar.ht",
      "params": {"reference_genome": "GRCh38", "export_tsv": true}
    },
    {
      "name": "interactome",
      "input": "/data/insider.bed.bgz",
      "output": "/out/interactome.ht",
      "params": {"reference_genome": "GRCh38"}
    }
  ]
}
```

Run:

```bash
hvantk mktable-batch --recipe examples/recipes/tables.example.json
```

YAML variant (examples/recipes/tables.example.yaml):

```yaml
---
tables:
  - name: clinvar
    input: /data/clinvar_2024.vcf.bgz
    output: /out/clinvar.ht
    params:
      reference_genome: GRCh38
      export_tsv: true
  - name: interactome
    input: /data/insider.bed.bgz
    output: /out/interactome.ht
    params:
      reference_genome: GRCh38
```

## 3) Build a single MatrixTable (MT)

- UCSC Cell Browser (TSV matrix + TSV metadata)

First, find and download the dataset you need:

```bash
# Search for datasets by tissue or keyword
hvantk download ucsc --list_datasets --search heart

# Download a dataset (use child path for collections)
hvantk download ucsc --dataset hoc/all-heart --output-dir data/ucsc
```

Then build the MatrixTable:

```bash
hvantk mkmatrix ucsc \
  --expression-matrix data/ucsc/hoc/all-heart/exprMatrix.tsv.gz \
  --metadata data/ucsc/hoc/all-heart/meta.tsv \
  --output-mt /out/ucsc.mt \
  --gene-column gene \
  --auto-convert-bgz \
  --overwrite
```

- Expression Atlas (TSV matrix + SDRF TSV)

```bash
hvantk mkmatrix expression-atlas \
  --expression-matrix /data/atlas/matrix.tsv \
  --sdrf /data/atlas/atlas.sdrf.tsv \
  --output-mt /out/atlas.mt \
  --gene-column "Gene ID" \
  --sample-id-column sample_id \
  --overwrite
```

If your expression matrix is a plain `.gz` file, add `--auto-convert-bgz`:

```bash
hvantk mkmatrix expression-atlas \
  --expression-matrix /data/atlas/matrix.tsv.gz \
  --sdrf /data/atlas/atlas.sdrf.tsv \
  --output-mt /out/atlas.mt \
  --auto-convert-bgz
```

- CPTAC (TSV/CSV expression + TSV/CSV metadata)

```bash
hvantk mkmatrix cptac \
  --expression /data/cptac/expression.tsv \
  --metadata /data/cptac/metadata.tsv \
  --output-mt /out/cptac.mt \
  --gene-id-col GeneID \
  --sample-id-col SampleID \
  --categorical-cols TumorType,Stage \
  --overwrite
```

## 4) Batch-create MatrixTables (MT) from a recipe

Example JSON recipe (save as examples/recipes/matrices.example.json):

```json
{
  "matrices": [
    {
      "name": "ucsc",
      "inputs": {
        "expression_matrix": "/data/ucsc/expr.tsv.bgz",
        "metadata": "/data/ucsc/meta.tsv"
      },
      "output": "/out/ucsc.mt",
      "params": {"gene_column": "gene", "overwrite": true}
    },
    {
      "name": "expression-atlas",
      "inputs": {
        "expression_matrix": "/data/atlas/matrix.tsv",
        "sdrf": "/data/atlas/atlas.sdrf.tsv"
      },
      "output": "/out/atlas.mt",
      "params": {"gene_column": "Gene ID", "sample_id_column": "sample_id"}
    }
  ]
}
```

Run:

```bash
hvantk mkmatrix-batch --recipe examples/recipes/matrices.example.json
```

YAML variant (examples/recipes/matrices.example.yaml):

```yaml
---
matrices:
  - name: ucsc
    inputs:
      expression_matrix: /data/ucsc/expr.tsv.bgz
      metadata: /data/ucsc/meta.tsv
    output: /out/ucsc.mt
    params:
      gene_column: gene
      overwrite: true
  - name: expression-atlas
    inputs:
      expression_matrix: /data/atlas/matrix.tsv
      sdrf: /data/atlas/atlas.sdrf.tsv
    output: /out/atlas.mt
    params:
      gene_column: "Gene ID"
      sample_id_column: sample_id
```

CPTAC JSON recipe (save as examples/recipes/cptac.example.json):

```json
{
  "matrices": [
    {
      "name": "cptac",
      "inputs": {
        "expression": "/data/cptac/expression.tsv",
        "metadata": "/data/cptac/metadata.tsv"
      },
      "output": "/out/cptac.mt",
      "params": {
        "gene_id_col": "GeneID",
        "gene_name_col": "Gene Name",
        "sample_id_col": "SampleID",
        "expression_col": "Expression",
        "categorical_cols": "TumorType,Stage",
        "overwrite": true
      }
    }
  ]
}
```

Run:

```bash
hvantk mkmatrix-batch --recipe examples/recipes/cptac.example.json
```

## 5) Ancestry Inference

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

Hail supports standard gzip (`.gz`) and uncompressed files but processes them single-threaded. Block gzip (BGZF) `.bgz` files enable parallel import and are strongly recommended for large datasets. hvantk provides two ways to convert:

### Automatic conversion with `--auto-convert-bgz`

Several CLI commands support the `--auto-convert-bgz` flag, which detects plain `.gz` files and converts them to BGZF before import:

```bash
hvantk mktable dbnsfp --raw-input data.gz --output-ht out.ht --auto-convert-bgz
hvantk mkmatrix ucsc -e expr.tsv.gz -m meta.tsv -o out.mt --auto-convert-bgz
hvantk mkmatrix expression-atlas -e matrix.tsv.gz -s atlas.sdrf.tsv -o out.mt --auto-convert-bgz
```

The converted `.bgz` file is written alongside the original (e.g., `data.gz` → `data.bgz`) and reused on subsequent runs.

### Standalone conversion with `utils convert-bgz`

For batch or one-off conversion:

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
