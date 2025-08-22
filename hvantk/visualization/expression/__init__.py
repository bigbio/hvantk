"""
Expression visualization namespace.

This subpackage groups plots for expression-related data across backends and data sources.
Current implementations:
- hvantk.visualization.expression.hail: Hail MatrixTable support

Additional backends (pandas/AnnData) may be added in future modules.
"""

# Intentionally avoid importing heavy backends here; import from submodules directly
