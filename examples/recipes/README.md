# Recipe Templates

JSON recipe templates for batch-building annotation tables and expression matrices.

For full documentation, see the [Recipes examples guide](https://bigbio.github.io/hvantk/examples/recipes/).

## Contents

| Recipe | Command | Description |
|--------|---------|-------------|
| `tables.example.json` | `hvantk mktable-batch` | ClinVar + INSIDER annotation tables |
| `matrices.example.json` | `hvantk mkmatrix-batch` | UCSC Cell Browser + Expression Atlas matrices |
| `cptac.example.json` | `hvantk mkmatrix-batch` | CPTAC protein expression matrices |

## Quick Start

```bash
# Build annotation tables
hvantk mktable-batch --recipe examples/recipes/tables.example.json

# Build expression matrices
hvantk mkmatrix-batch --recipe examples/recipes/matrices.example.json

# Build CPTAC protein matrices
hvantk mkmatrix-batch --recipe examples/recipes/cptac.example.json
```

## Customizing

Copy a template and edit paths and parameters:

```bash
cp examples/recipes/tables.example.json my_recipe.json
# Edit input/output paths in my_recipe.json
hvantk mktable-batch --recipe my_recipe.json
```
