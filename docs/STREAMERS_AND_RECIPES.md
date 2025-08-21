# Streamers and Recipes

This document defines a small, testable plugin interface for streamers and a JSON recipe format (YAML optional) to compose pipelines.

## Streamer interface (Hail)

- Base classes: hvantk/data/data_streamer.py
  - DataStreamer: setup(), stream(), process_chunk(), teardown()
  - HailDataStreamer: wraps Hail init/stop
- Annotation streamers: hvantk/annotation/annotation_streamer.py provides examples like:
  - VariantPredictionScoreStreamer
  - GeneExpressionStreamer
  - GeneConstraintStreamer
  - PopulationFrequencyStreamer

Guidelines:
- Stateless: no hidden global state; take all config via __init__
- Small: aim for 50–150 lines, single responsibility
- Chunk-aware: stream() yields small HTs/MTs or process_chunk() handles them

## Recipe format (JSON by default)

A recipe wires inputs, selectors, and streamers to produce outputs reproducibly.

Fields:
- name: string
- inputs:
  - type: table|matrixtable
  - path: local path or catalog id
  - selectors: optional dict (genes, regions, tissues, cell_types, timepoints, cohorts)
- streamers: list of steps; each has:
  - name: class name of streamer
  - params: key-value pairs passed to __init__
- output:
  - path: destination (HT/MT)
  - export: optional list [tsv, parquet]

Example (JSON):

```json
{
  "name": "heart-rare-variants",
  "inputs": {
    "type": "table",
    "path": "catalog:clinvar.ht@2025-05-01",
    "selectors": {
      "genes_panel": "data/panels/chd_genes.txt"
    }
  },
  "streamers": [
    { "name": "PopulationFrequencyStreamer", "params": {} },
    { "name": "GeneConstraintStreamer", "params": {} },
    { "name": "GeneExpressionStreamer", "params": { "tissue_focus": "heart" } }
  ],
  "output": {
    "path": "data/recipes/heart_rare_variants.ht",
    "export": ["tsv"]
  }
}
```

(YAML is also supported if PyYAML is installed; see examples/recipes/heart_rare_variants.yaml for a YAML variant.)

## Minimal Python runner

Until a CLI exists, you can prototype a recipe runner in a notebook or script:

- Load the input HT/MT (from path or by resolving a catalog id; see DATA_CATALOG.md)
- Instantiate streamers by name and params
- Apply each streamer in order using process_chunk()
- Write/checkpoint the result, then export if requested

See examples for a tiny JSON recipe (YAML optional) and how to bind it to available streamers.
