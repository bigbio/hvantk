# AlphaGenome Variant Prediction Streamer — Design Spec

**Date:** 2026-03-31
**Branch:** `variantpredictions`
**Status:** Approved

## Overview

Add an `AlphaGenomeStreamer` to hvantk that, given a set of variant positions, calls the AlphaGenome API and produces per-modality Hail Tables with the full multimodal predictions (expression, splicing, chromatin, contact maps). The streamer handles authentication, adaptive interval grouping, rate limiting, checkpoint-based resumption, and assembly of final outputs.

## Motivation

AlphaGenome (DeepMind, 2026) is a genomic foundation model that predicts regulatory effects of genetic variants at single base-pair resolution across multiple modalities. Integrating it into hvantk enables researchers to enrich variant annotations with predicted functional impact — complementing existing static annotation sources (ClinVar, gnomAD, dbNSFP) with model-based predictions.

## Config File

A YAML config file controls API access, ontology terms, and interval strategy:

```yaml
api:
  key: "AG-xxxx"                    # or null -> falls back to ALPHAGENOME_API_KEY env var
  max_retries: 3
  retry_backoff: 2.0                # exponential backoff base (seconds)
  request_timeout: 120              # seconds per API call

ontology:
  terms:
    - "UBERON:0001157"              # liver
    - "UBERON:0000955"              # brain
    - "UBERON:0002048"              # lung
  output_types:
    - RNA_SEQ
    - CHROMATIN
    - SPLICING
    - CONTACT_MAP

intervals:
  default_size: 1048576             # 1Mbp fixed-size fallback
  adaptive: true                    # group nearby variants into shared intervals
  adaptive_max_size: 1048576        # max interval size when adaptive
  density_window: 50000             # bp window to check for neighboring variants
```

**Auth resolution order:** config `api.key` > `ALPHAGENOME_API_KEY` env var > error with message.

## Interval Strategy

### Adaptive Mode (default)

1. Sort all variants by genomic position (chrom, pos).
2. Scan linearly — group consecutive variants within `density_window` (default 50kb) of each other.
3. For each group, compute a bounding interval centered on the group midpoint, expanded to cover all variants plus padding.
4. Cap at `adaptive_max_size` (1Mbp — AlphaGenome's limit). If a group exceeds this, split into sub-groups.
5. Singleton variants (no neighbors within the window) get `default_size` centered on their position.

### Fixed Mode (`adaptive: false`)

Each variant gets a `default_size` window centered on its position. No grouping — one API call per variant.

### Benefit

If 200 of 5000 variants cluster in a gene-dense region, they share a single API call instead of 200 separate ones. For a typical GWAS-like variant set this can reduce API calls significantly.

### Implementation

A pure function `compute_intervals(variants, config) -> List[Tuple[Interval, List[Variant]]]` that returns interval-to-variant mappings. Testable independently of the API.

## Streamer Architecture

`AlphaGenomeStreamer` extends `HailDataStreamer`.

### Input

- Hail Table keyed by `(locus, alleles)`, or
- TSV with `chrom, pos, ref, alt` columns (auto-imported into Hail Table on `setup()`).

### Lifecycle

#### `setup()`

1. Initialize Hail.
2. Load and validate config YAML.
3. Authenticate — create `dna_client` via `dna_client.create(api_key)`.
4. Load input variants (Hail Table or TSV).
5. Load checkpoint state if resuming — a JSON file tracking completed intervals.
6. Compute interval-to-variant mappings via `compute_intervals()`.
7. Filter out already-completed intervals from checkpoint.

#### `stream()`

- Iterate over pending intervals in batches of `chunk_size`.
- For each interval + its variants, call `model.predict_variant()` per variant.
- Collect raw outputs into per-modality dicts.
- After each batch: write intermediate results to checkpoint dir, update checkpoint JSON.
- Yield a dict of `{modality_name: hl.Table}` per batch (intermediate, used by `StreamProcessor` pipeline or discarded if running standalone).

#### `teardown()`

- Union all batch checkpoints into final per-modality Hail Tables.
- Write to output paths (e.g., `output_dir/rna_seq.ht`, `output_dir/chromatin.ht`).
- Each table keyed by `(locus, alleles)` with ontology term as a column/nested struct.
- Log summary stats (variants processed, API calls made, failures).

### Checkpoint Directory Structure

```
output_dir/
├── _checkpoints/
│   ├── state.json          # {"completed_intervals": [...], "failed_variants": [...]}
│   ├── batch_000.json      # raw API responses for batch 0
│   ├── batch_001.json
│   └── ...
├── rna_seq.ht/
├── chromatin.ht/
├── splicing.ht/
└── contact_map.ht/
```

## Rate Limiting & Error Handling

### Adaptive Throttling

AlphaGenome doesn't publish fixed rate limits (they vary by demand), so the streamer uses adaptive throttling:

- **Base delay:** Configurable minimum wait between API calls (default 0.5s).
- **Backoff on 429/5xx:** Exponential backoff with jitter — `retry_backoff^attempt * (1 + random(0, 0.5))`, up to `max_retries`.
- **Cooldown:** If 3 consecutive requests get rate-limited, pause for 60s before continuing.
- **Progress logging:** Log every N intervals completed, estimated time remaining based on average call duration.

### Error Classification

| Error Type | Behavior |
|------------|----------|
| 429 Too Many Requests | Backoff + retry |
| 5xx Server Error | Backoff + retry |
| 4xx Client Error (not 429) | Log warning, skip variant, record in `failed_variants` |
| Timeout | Retry up to `max_retries`, then skip + record |
| Network Error | Retry up to `max_retries`, then skip + record |
| Invalid variant (no ref/alt at position) | Log warning, skip, record |

### Resumption Flow

1. On start, check for `_checkpoints/state.json`.
2. If exists, load completed intervals — skip them, log "Resuming from batch N".
3. Failed variants from previous run are re-attempted once.
4. User can force a clean start via `--no-resume` CLI flag.

### Final Report

After completion, log a summary: total variants, successful, failed (with reasons), API calls made, total runtime.

## CLI & Table Builder Integration

### CLI Command

```bash
hvantk mktable alphagenome \
  --input variants.ht \
  --output-dir /out/alphagenome/ \
  --config alphagenome_config.yaml \
  --no-resume \
  --overwrite
```

### Builder Function

In `hvantk/tables/table_builders.py`:

```python
def create_alphagenome_tb(
    input_path: str,
    output_path: str,
    config_path: str,
    no_resume: bool = False,
    overwrite: bool = False,
) -> Dict[str, hl.Table]:
    """Build per-modality Hail Tables from AlphaGenome predictions."""
```

Returns `Dict[str, hl.Table]` — a deviation from the standard single-table builder signature, necessary for multi-table output.

### Registry

In `hvantk/tables/registry.py`:

```python
TABLE_BUILDERS["alphagenome"] = create_table_adapter(
    "hvantk.tables.table_builders", "create_alphagenome_tb"
)
```

### Recipe Support

Output is a directory path:

```json
{
  "name": "alphagenome",
  "input": "/data/variants.ht",
  "output": "/out/alphagenome/",
  "params": {"config_path": "alphagenome_config.yaml"}
}
```

## Files & Module Layout

### New Files

- `hvantk/data/alphagenome_streamer.py` — `AlphaGenomeStreamer` class + `compute_intervals()` function
- `hvantk/tests/test_alphagenome_streamer.py` — unit tests with mocked API
- `hvantk/tests/testdata/alphagenome_config.yaml` — test config fixture

### Modified Files

- `hvantk/tables/table_builders.py` — add `create_alphagenome_tb()`
- `hvantk/tables/registry.py` — register `"alphagenome"` builder
- `hvantk/commands/make_table_cli.py` — add `mktable_alphagenome` command
- `hvantk/core/constants.py` — add AlphaGenome default constants
- `hvantk/data/__init__.py` — export `AlphaGenomeStreamer`

### Dependencies

- `alphagenome` (PyPI) — **optional dependency**. Builder raises `ImportError` with install instructions if missing.
- No new core dependencies.

## Testing Strategy

- Unit tests with mocked `dna_client` — test streamer lifecycle, interval computation, checkpoint logic, error handling.
- `compute_intervals()` tested as a pure function with known variant distributions.
- Integration tests marked `@pytest.mark.network` for real API calls (small variant set).
- Config validation tests (missing key, invalid ontology terms, bad interval params).
