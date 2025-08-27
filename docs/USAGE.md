# hvantk Usage Guide

This guide shows practical, copy-pasteable examples to build Hail Tables (HT) and MatrixTables (MT) from explicit raw files and from recipes (JSON/YAML).

If you haven’t installed hvantk yet, see the main README for install steps.

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

- Ensembl gene annotations (Biomart TSV keyed by gene_id)

```bash
hvantk mktable ensembl-gene \
  --raw-input /data/biomart.tsv.bgz \
  --output-ht /out/ensembl.ht \
  --no-canonical
```

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

```bash
hvantk mkmatrix ucsc \
  --expression-matrix /data/ucsc/expr.tsv.bgz \
  --metadata /data/ucsc/meta.tsv \
  --output-mt /out/ucsc.mt \
  --gene-column gene \
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

## Tips & troubleshooting

- Use `--overwrite` to replace an existing output. Without it, builders abort if the output exists.
- For JSON vs YAML: JSON works out of the box; YAML recipes require `PyYAML` installed.
- For UCSC, gene labels may be pipe-delimited (e.g., A|B); `--split-gene-field` defaults to true.
- MatrixTables typically store sample/cell metadata under `mt.col_key` and cols metadata; inspect with `mt.describe()` in Python or logs from CLI.

## See also

- docs/DATA_CATALOG.md – versioning, provenance, and hosting strategy
- docs/DEVELOPING.md – dev workflow and builder contracts
- examples/recipes/ – ready-to-edit recipe templates
