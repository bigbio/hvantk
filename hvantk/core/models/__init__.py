# --- Phase A artifact exports ---
from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.expression_matrix import ExpressionMatrix
from hvantk.core.models.gene_set import GeneSet
from hvantk.core.models.provenance import Provenance
from hvantk.core.models.build_context import BuildContext
from hvantk.core.models._expr import col, agg_count, agg_mean, agg_sum

# Prevent underscore submodule names from leaking into the public namespace.
# Python injects them as package attributes when importing from them;
# deleting here keeps dir(hvantk.core.models) clean.
try:
    del _expr  # type: ignore[name-defined]
except NameError:
    pass
try:
    del _compile  # type: ignore[name-defined]
except NameError:
    pass

__all__ = [
    "AnnotationTable",
    "ExpressionMatrix",
    "GeneSet",
    "Provenance",
    "BuildContext",
    "col",
    "agg_count",
    "agg_mean",
    "agg_sum",
]
