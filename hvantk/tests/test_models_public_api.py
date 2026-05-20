"""The public API surface of hvantk.core.models is what plugin authors
and algorithms see. Pinning it as a test catches accidental removals.
"""
from __future__ import annotations


def test_public_exports_available():
    from hvantk.core.models import (
        AnnotationTable,
        ExpressionMatrix,
        GeneSet,
        Provenance,
        BuildContext,
        col,
        agg_count,
        agg_mean,
        agg_sum,
    )
    assert AnnotationTable is not None
    assert ExpressionMatrix is not None
    assert GeneSet is not None
    assert Provenance is not None
    assert BuildContext is not None
    assert callable(col)
    assert callable(agg_count)
    assert callable(agg_mean)
    assert callable(agg_sum)


def test_underscore_modules_not_in_dir():
    """Underscore-prefixed modules are implementation detail — don't re-export them."""
    import hvantk.core.models as m

    assert "_expr" not in dir(m)
    assert "_compile" not in dir(m)
