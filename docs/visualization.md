# Visualization package refactor

This document describes the structure of the `hvantk.visualization` package and migration notes.

## Goals
- Scalable layout for rich, multi-omics visualization
- Clear separation by data domain and backend
- Backward compatibility for existing imports

## Package layout
- `hvantk.visualization.base` — shared helpers (style, save, colors, figure labels)
- `hvantk.visualization.expression` — expression-focused plots
  - `hvantk.visualization.expression.hail` — Hail MatrixTable implementations
- (planned) `hvantk.visualization.variants` — variant/VCF/MatrixTable plots
- (planned) `hvantk.visualization.multiomics` — cross-omics dashboards/composites
- (planned) `hvantk.visualization.interactive` — Plotly/Bokeh/Altair wrappers
- (planned) `hvantk.visualization.qc` — generic QC/diagnostics

## Public API and imports
High-level functions are re-exported from `hvantk.visualization` using lazy loading to avoid heavy imports at module import time.

Current public function:
- `visualize_expression_distribution(mt, n_bins=50, log_scale=True, title=..., expr_field='x')`
  - Import paths (both supported):
    - `from hvantk.visualization import visualize_expression_distribution`
    - `from hvantk.visualization.expression.hail import visualize_expression_distribution`

## Deprecations
- `hvantk.visualization.hail_expression` is deprecated and will be removed in a future release.
  - A shim remains in place and will emit a `DeprecationWarning` while forwarding imports.
  - Update imports to one of the new paths above.

## Notes
- The package-level re-export is a thin wrapper that imports the backend on first use, keeping `import hvantk.visualization` lightweight.
- Additional backends (pandas/AnnData) can be added as sibling modules under `hvantk.visualization.expression` without breaking the public API.

