# Streamers and Recipes

This document describes the streamer interface and recipe format for hvantk data processing pipelines.

## Overview

Streamers are stateless, composable data transformers that follow the `Streamer` protocol:
- **Transform**: Apply operations to Hail data structures (filtering, annotation, aggregation, joining)
- **Validate**: Check that input data matches expected schema
- **Metadata**: Provide information about the streamer's capabilities and requirements

Streamers focus purely on transformation logic. I/O operations (reading from sources and writing to destinations) 
are handled by the pipeline orchestration layer, which allows streamers to remain simple, testable, and composable.

Recipes are JSON or YAML configurations that chain streamers together to answer biological questions.

## Streamer Interface

A streamer should implement the `Streamer` protocol defined in `hvantk/core/protocols.py`:

```python
from typing import Union, Dict, Any
import hail as hl

class Streamer:
    def transform(
        self, 
        input_data: Union[hl.Table, hl.MatrixTable], 
        **params: Any
    ) -> Union[hl.Table, hl.MatrixTable]:
        """
        Transform input data and return output.
        
        Parameters
        ----------
        input_data : Union[hl.Table, hl.MatrixTable]
            Input Hail data structure
        **params : Any
            Transformation parameters (e.g., filter thresholds, join tables, etc.)
        
        Returns
        -------
        Union[hl.Table, hl.MatrixTable]
            Transformed Hail data structure
        """
        pass
    
    def validate_input(
        self, 
        input_data: Union[hl.Table, hl.MatrixTable]
    ) -> bool:
        """
        Validate that input schema matches expectations.
        
        Parameters
        ----------
        input_data : Union[hl.Table, hl.MatrixTable]
            The Hail data structure to validate
        
        Returns
        -------
        bool
            True if input schema is valid, False otherwise
        """
        pass
    
    def get_metadata(self) -> Dict[str, Any]:
        """
        Return metadata about this streamer.
        
        Returns
        -------
        Dict[str, Any]
            Metadata dictionary containing:
            - type: str - 'filter', 'join', 'aggregate', 'annotate', etc.
            - input_type: str - Expected input type
            - output_type: str - Output type (optional)
            - description: str - Human-readable description (optional)
        """
        pass
```

**Note**: Streamers focus on transformation logic only. Reading from sources and writing to destinations 
are handled separately by pipeline orchestration (see recipe execution below). This separation allows 
streamers to be pure, stateless transformers that are easily testable and composable.

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

## Built-in Pipeline Components

The pipeline supports various components for different stages of data processing:

### I/O Components (Readers/Writers)
These handle data loading and exporting, separate from the Streamer protocol:
- `table_reader`: Load Hail Table
- `matrix_table_reader`: Load Hail MatrixTable
- `vcf_reader`: Load VCF file
- `table_writer`: Write Hail Table
- `matrix_table_writer`: Write Hail MatrixTable
- `vcf_writer`: Export to VCF
- `tsv_writer`: Export to TSV

### Streamers (Transform Components)
These implement the `Streamer` protocol for data transformation:
- `frequency_filter`: Filter by allele frequency
- `quality_filter`: Filter by quality metrics
- `expression_annotator`: Add expression data
- `consequence_annotator`: Add variant consequences

## Creating Custom Streamers

To create a custom streamer, implement the `Streamer` protocol:

1. Implement all protocol methods: `transform()`, `validate_input()`, and `get_metadata()`
2. Register it in the streamer registry
3. Add tests with fixtures
4. Document parameters

Example:

```python
from typing import Union, Dict, Any
import hail as hl

class MyCustomStreamer:
    """Custom streamer that filters variants by a custom criterion."""
    
    def transform(
        self, 
        input_data: Union[hl.Table, hl.MatrixTable], 
        **params: Any
    ) -> Union[hl.Table, hl.MatrixTable]:
        """Apply custom filtering logic."""
        threshold = params.get('threshold', 0.5)
        filtered = input_data.filter(input_data.custom_field >= threshold)
        return filtered
    
    def validate_input(
        self, 
        input_data: Union[hl.Table, hl.MatrixTable]
    ) -> bool:
        """Validate that input has required field."""
        return 'custom_field' in input_data.row
    
    def get_metadata(self) -> Dict[str, Any]:
        """Return streamer metadata."""
        return {
            'type': 'filter',
            'input_type': 'variant_table',
            'description': 'Filters variants by custom field threshold'
        }

# Register the streamer
from hvantk.streamers import register_streamer
register_streamer('my_custom_filter', MyCustomStreamer)
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
