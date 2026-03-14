# Recipe Templates

hvantk supports batch processing of annotation tables and expression matrices using JSON recipe files. Recipes let you define multiple build tasks in a single file and execute them all at once.

## How Recipes Work

A recipe is a JSON file that specifies a list of tables or matrices to build, along with their input paths, output paths, and optional parameters. You run recipes with the `mktable-batch` or `mkmatrix-batch` CLI commands:

```bash
# Build annotation tables from a recipe
hvantk mktable-batch --recipe path/to/recipe.json

# Build expression matrices from a recipe
hvantk mkmatrix-batch --recipe path/to/recipe.json
```

All example recipe files are available in the repository at [`examples/recipes/`](https://github.com/bigbio/hvantk/tree/main/examples/recipes/).

## Table Recipes

**[`tables.example.json`](https://github.com/bigbio/hvantk/tree/main/examples/recipes/tables.example.json)** defines annotation tables (Hail Tables) to build in batch.

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

**Run it:**

```bash
hvantk mktable-batch --recipe tables.example.json
```

This demonstrates:
- ClinVar annotation table creation from a VCF file
- INSIDER interactome table creation from a BED file
- Reference genome specification (`GRCh38`)
- Optional TSV export (`export_tsv: true`)

## Matrix Recipes

**[`matrices.example.json`](https://github.com/bigbio/hvantk/tree/main/examples/recipes/matrices.example.json)** defines expression matrices (Hail MatrixTables) to build in batch.

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

**Run it:**

```bash
hvantk mkmatrix-batch --recipe matrices.example.json
```

This demonstrates:
- UCSC Cell Browser data conversion (expression matrix + metadata)
- Expression Atlas data processing (expression matrix + SDRF)
- Multiple matrix creation in a single run

## CPTAC Recipes

**[`cptac.example.json`](https://github.com/bigbio/hvantk/tree/main/examples/recipes/cptac.example.json)** defines CPTAC protein expression matrices.

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
        "categorical_cols": ["TumorType", "Stage"],
        "numeric_cols": ["Age", "Expression"],
        "overwrite": true
      }
    }
  ]
}
```

**Run it:**

```bash
hvantk mkmatrix-batch --recipe cptac.example.json
```

This demonstrates:
- Protein expression matrix creation from CPTAC data
- Categorical metadata handling (tumor type, stage)
- Custom column mapping for gene and sample identifiers
- Numeric column specification

## Creating Your Own Recipes

1. Copy an example recipe as a starting point:
   ```bash
   cp examples/recipes/tables.example.json my_recipe.json
   ```

2. Edit paths and parameters to match your data:
   ```json
   {
     "tables": [
       {
         "name": "clinvar",
         "input": "/path/to/your/clinvar.vcf.bgz",
         "output": "/path/to/your/output.ht",
         "params": {"reference_genome": "GRCh38"}
       }
     ]
   }
   ```

3. Run your recipe:
   ```bash
   hvantk mktable-batch --recipe my_recipe.json
   ```

## Recipe Format Reference

### Table Recipe Fields

| Field | Required | Description |
|-------|----------|-------------|
| `name` | Yes | Builder name (e.g., `clinvar`, `interactome`, `dbnsfp`) |
| `input` | Yes | Path to the input file |
| `output` | Yes | Path for the output Hail Table (`.ht`) |
| `params` | No | Builder-specific parameters (varies by source) |

### Matrix Recipe Fields

| Field | Required | Description |
|-------|----------|-------------|
| `name` | Yes | Builder name (e.g., `ucsc`, `expression-atlas`, `cptac`) |
| `inputs` | Yes | Object with named input paths |
| `output` | Yes | Path for the output MatrixTable (`.mt`) |
| `params` | No | Builder-specific parameters (varies by source) |

## Further Reading

For full recipe documentation including all supported builders and their parameters, see the [Usage Guide](../guide/usage.md).
