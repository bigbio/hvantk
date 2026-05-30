# --- Phase A artifact exports ---
from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.expression_matrix import ExpressionMatrix
from hvantk.core.models.gene_set import GeneSet
from hvantk.core.models.variant_matrix import VariantMatrix
from hvantk.core.models.provenance import Provenance
from hvantk.core.models.build_context import BuildContext
from hvantk.core.models._expr import col, agg_count, agg_mean, agg_sum

# Prevent underscore submodule names from leaking into the public namespace.
# Python injects them as package attributes when importing from them;
# popping from globals() keeps dir(hvantk.core.models) clean. We use
# globals().pop() instead of `del _expr` so flake8 doesn't flag the
# names as undefined (they only exist as side-effects of the imports above).
globals().pop("_expr", None)
globals().pop("_compile", None)

__all__ = [
    "AnnotationTable",
    "ExpressionMatrix",
    "GeneSet",
    "VariantMatrix",
    "Provenance",
    "BuildContext",
    "col",
    "agg_count",
    "agg_mean",
    "agg_sum",
]
