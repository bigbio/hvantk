"""Typed artifacts and the build-time helpers that stamp them.

Every name below is resolved lazily (PEP 562). Importing them eagerly made the
*package* import pull in pandas (via ``annotation_table``) and anndata (via
``expression_matrix``) -- roughly 1.2 s -- and Python imports a package before
any of its submodules, so even ``from hvantk.core.models.provenance import
Provenance`` (a plain dataclass) paid that cost. The plugin loader imports
exactly that, which is why every ``hvantk plugins``/``catalog``/``reprocess``
invocation used to wait on anndata before printing anything.

The public surface is unchanged: ``from hvantk.core.models import
AnnotationTable`` still works, it just imports on first use.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - import-time typing only
    from hvantk.core.models.annotation_table import AnnotationTable
    from hvantk.core.models.build_context import BuildContext
    from hvantk.core.models.expression_matrix import ExpressionMatrix
    from hvantk.core.models.gene_set import GeneSet
    from hvantk.core.models.provenance import Provenance
    from hvantk.core.models.variant_matrix import VariantMatrix
    from hvantk.core.models._expr import agg_count, agg_mean, agg_sum, col

# public name -> defining submodule
_EXPORTS = {
    "AnnotationTable": "hvantk.core.models.annotation_table",
    "ExpressionMatrix": "hvantk.core.models.expression_matrix",
    "GeneSet": "hvantk.core.models.gene_set",
    "VariantMatrix": "hvantk.core.models.variant_matrix",
    "Provenance": "hvantk.core.models.provenance",
    "BuildContext": "hvantk.core.models.build_context",
    "col": "hvantk.core.models._expr",
    "agg_count": "hvantk.core.models._expr",
    "agg_mean": "hvantk.core.models._expr",
    "agg_sum": "hvantk.core.models._expr",
}

__all__ = [
    "AnnotationTable",
    "BuildContext",
    "ExpressionMatrix",
    "GeneSet",
    "Provenance",
    "VariantMatrix",
    "agg_count",
    "agg_mean",
    "agg_sum",
    "col",
]

if set(__all__) != set(_EXPORTS):  # not an assert: python -O strips those
    raise RuntimeError(
        "hvantk.core.models: __all__ and _EXPORTS disagree -- "
        f"{sorted(set(__all__) ^ set(_EXPORTS))}"
    )


def __getattr__(name: str):
    module = _EXPORTS.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib

    value = getattr(importlib.import_module(module), name)
    globals()[name] = value  # cache, so the import cost is paid once
    return value


def __dir__():
    # Union of the lazy exports and whatever is already bound, minus private
    # names. The underscore filter is what the old `globals().pop("_expr")`
    # achieved; returning only __all__ would additionally hide the real
    # submodules (annotation_table, provenance, ...) from dir().
    return sorted(
        {*__all__, *(n for n in globals() if not n.startswith("_"))} - {"TYPE_CHECKING"}
    )
