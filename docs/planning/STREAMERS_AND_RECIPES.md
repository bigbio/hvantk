# Streamers and Recipes

This document describes the streamer interface and recipe format for hvantk data processing pipelines.

## Overview

Streamers are stateless, composable data transformers that follow a simple contract:
- **Read**: Load data from a source (Hail Table, MatrixTable, or file)
- **Transform**: Apply operations (filtering, annotation, aggregation)
- **Write**: Output results to a destination

Recipes are JSON or YAML configurations that chain streamers together to answer biological questions.

## Streamer Interface

A streamer should implement:

```python
class Streamer:
    def read(self, source: str) -> Union[hl.Table, hl.MatrixTable]:
        """Load data from source"""
        pass
    
    def transform(self, data: Union[hl.Table, hl.MatrixTable]) -> Union[hl.Table, hl.MatrixTable]:
        """Apply transformations"""
        pass
    
    def write(self, data: Union[hl.Table, hl.MatrixTable], destination: str) -> None:
        """Write data to destination"""
        pass
```

### Design Principles

1. **Stateless**: No internal state between operations
2. **Composable**: Output of one streamer is valid input for another
3. **Testable**: Can be tested on tiny fixtures
4. **Parameterized**: Accept configuration via arguments

## Recipe Format

Recipes define data processing pipelines using JSON or YAML.

### JSON Recipe Example

```json
{
  "name": "rare_variant_analysis",
  "description": "Identify rare variants in tissue-specific genes",
  "steps": [
    {
      "name": "load_variants",
      "type": "table_reader",
      "params": {
        "input": "/data/variants.ht"
      }
    },
    {
      "name": "filter_rare",
      "type": "frequency_filter",
      "params": {
        "max_af": 0.01,
        "population": "global"
      }
    },
    {
      "name": "annotate_expression",
      "type": "expression_annotator",
      "params": {
        "expression_mt": "/data/heart_expression.mt",
        "tissue": "heart",
        "min_expression": 5.0
      }
    },
    {
      "name": "export_results",
      "type": "table_writer",
      "params": {
        "output": "/results/rare_heart_variants.ht",
        "format": "hail_table"
      }
    }
  ]
}
```

### YAML Recipe Example

```yaml
---
name: rare_variant_analysis
description: Identify rare variants in tissue-specific genes
steps:
  - name: load_variants
    type: table_reader
    params:
      input: /data/variants.ht
  
  - name: filter_rare
    type: frequency_filter
    params:
      max_af: 0.01
      population: global
  
  - name: annotate_expression
    type: expression_annotator
    params:
      expression_mt: /data/heart_expression.mt
      tissue: heart
      min_expression: 5.0
  
  - name: export_results
    type: table_writer
    params:
      output: /results/rare_heart_variants.ht
      format: hail_table
```

## Recipe Execution

Recipes are executed sequentially, with each step's output becoming the next step's input:

```bash
# Execute a recipe
hvantk run-recipe --recipe /path/to/recipe.json

# Or in Python
from hvantk.recipes import execute_recipe
execute_recipe('recipe.json')
```

## Built-in Streamers

### Readers
- `table_reader`: Load Hail Table
- `matrix_table_reader`: Load Hail MatrixTable
- `vcf_reader`: Load VCF file

### Transformers
- `frequency_filter`: Filter by allele frequency
- `quality_filter`: Filter by quality metrics
- `expression_annotator`: Add expression data
- `consequence_annotator`: Add variant consequences

### Writers
- `table_writer`: Write Hail Table
- `matrix_table_writer`: Write Hail MatrixTable
- `vcf_writer`: Export to VCF
- `tsv_writer`: Export to TSV

## Creating Custom Streamers

To create a custom streamer:

1. Implement the streamer interface
2. Register it in the streamer registry
3. Add tests with fixtures
4. Document parameters

Example:

```python
from hvantk.streamers import BaseStreamer

class MyCustomStreamer(BaseStreamer):
    def transform(self, data, **kwargs):
        # Apply custom logic
        filtered = data.filter(...)
        return filtered

# Register
from hvantk.streamers import register_streamer
register_streamer('my_custom', MyCustomStreamer)
```

## Testing Recipes

Test recipes on small fixtures:

```bash
# Use test data
hvantk run-recipe --recipe recipe.json --test-mode

# Validate without executing
hvantk validate-recipe --recipe recipe.json
```

## Best Practices

1. **Keep streamers focused**: Each streamer should do one thing well
2. **Use parameter caching**: Cache results by parameter hash
3. **Support selectors**: Allow filtering to genes/regions/tissues before building
4. **Document clearly**: Include parameter descriptions and examples
5. **Version recipes**: Track recipe versions with data versions

## See Also

- [Data Catalog](DATA_CATALOG.md) - Dataset registry and versioning
- [Developer Guide](DEVELOPING.md) - Development workflow
- [Usage Guide](../library/usage.md) - End-user documentation
- [Examples](../../examples/recipes/) - Ready-to-use recipe templates
