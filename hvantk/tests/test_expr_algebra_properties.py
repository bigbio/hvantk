"""Phase O: property-based tests for the Expr algebra.

These tests assert algebraic laws that must hold for any expression on the
pandas backend.  They catch subtle compile_to_pandas divergences that the
per-operator parity tests in test_expr_compile.py would miss.

Hypothesis is used for random data generation with ``max_examples=20`` and
``deadline=None`` to keep runtime predictable (typically <30 s total).

All strategies avoid generating null values so pandas NaN propagation does
not interfere with boolean-algebra identities that assume a two-valued logic.
"""
from __future__ import annotations

from datetime import datetime, timezone

import pandas as pd
import pytest
from hypothesis import given, settings, strategies as st

from hvantk.core.models._expr import col, Expr  # noqa: F401
from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.provenance import Provenance


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _prov() -> Provenance:
    return Provenance(
        plugin="t",
        dataset="t:rows",
        plugin_version="0",
        source_fingerprint="sha256:x",
        schema_id="t-rows-v1",
        build_timestamp=datetime(2026, 5, 21, tzinfo=timezone.utc),
        builder_commit=None,
    )


def _make_ann(rows: list[dict]) -> AnnotationTable:
    return AnnotationTable.from_pandas(pd.DataFrame(rows), provenance=_prov())


def _sort_rows(rows: list[dict]) -> list[dict]:
    """Return rows in a stable order so set-equality checks work."""
    return sorted(rows, key=lambda r: tuple(sorted((k, str(v)) for k, v in r.items())))


# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

# Column ``x``: small ints (no nulls)
_int_strat = st.integers(min_value=-50, max_value=50)

# Column ``s``: short strings from a fixed alphabet (no nulls)
_str_val_strat = st.sampled_from(["a", "b", "c", "d"])


@st.composite
def _small_ann(draw) -> AnnotationTable:
    """Generate a small AnnotationTable with columns x (int), y (int), s (str)."""
    n = draw(st.integers(min_value=1, max_value=20))
    xs = draw(st.lists(_int_strat, min_size=n, max_size=n))
    ys = draw(st.lists(st.integers(min_value=0, max_value=10), min_size=n, max_size=n))
    ss = draw(st.lists(_str_val_strat, min_size=n, max_size=n))
    rows = [{"x": x, "y": y, "s": s} for x, y, s in zip(xs, ys, ss)]
    return _make_ann(rows)


# Predicate factories (callables so each draw creates a fresh Expr node)
_PRED_FACTORIES = [
    lambda: col("x") > 0,
    lambda: col("x") < 50,
    lambda: col("x") == 0,
    lambda: col("y") >= 5,
    lambda: col("s") == "a",
    lambda: col("s").isin(["a", "b"]),
    lambda: col("x") != 0,
    lambda: col("y") < 8,
]

_pred_strat = st.sampled_from(_PRED_FACTORIES)

_SETTINGS = dict(max_examples=20, deadline=None)


# ---------------------------------------------------------------------------
# 1. Filter composition law: filter(a).filter(b) ≡ filter(a & b)
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann(), pred_a=_pred_strat, pred_b=_pred_strat)
def test_filter_composition_equiv_to_and(ann, pred_a, pred_b):
    """filter(a).filter(b) ≡ filter(a & b)"""
    a, b = pred_a(), pred_b()
    chained = _sort_rows(ann.filter(a).filter(b).collect())
    combined = _sort_rows(ann.filter(a & b).collect())
    assert chained == combined


# ---------------------------------------------------------------------------
# 2. Double-negation: filter(~~a) ≡ filter(a)
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann(), pred=_pred_strat)
def test_filter_double_negation(ann, pred):
    """filter(~~a).collect() ≡ filter(a).collect()"""
    p = pred()
    assert _sort_rows(ann.filter(~~p).collect()) == _sort_rows(ann.filter(p).collect())


# ---------------------------------------------------------------------------
# 3. And absorption: filter(a & a) ≡ filter(a)
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann(), pred=_pred_strat)
def test_filter_and_absorption(ann, pred):
    """filter(a & a) ≡ filter(a)"""
    p = pred()
    assert _sort_rows(ann.filter(p & p).collect()) == _sort_rows(ann.filter(p).collect())


# ---------------------------------------------------------------------------
# 4. Or absorption: filter(a | a) ≡ filter(a)
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann(), pred=_pred_strat)
def test_filter_or_absorption(ann, pred):
    """filter(a | a) ≡ filter(a)"""
    p = pred()
    assert _sort_rows(ann.filter(p | p).collect()) == _sort_rows(ann.filter(p).collect())


# ---------------------------------------------------------------------------
# 5. And commutativity: filter(a & b) ≡ filter(b & a)
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann(), pred_a=_pred_strat, pred_b=_pred_strat)
def test_and_commutative(ann, pred_a, pred_b):
    """filter(a & b) ≡ filter(b & a)  (as multisets)"""
    a, b = pred_a(), pred_b()
    left = _sort_rows(ann.filter(a & b).collect())
    right = _sort_rows(ann.filter(b & a).collect())
    assert left == right


# ---------------------------------------------------------------------------
# 6. Or commutativity: filter(a | b) ≡ filter(b | a)
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann(), pred_a=_pred_strat, pred_b=_pred_strat)
def test_or_commutative(ann, pred_a, pred_b):
    """filter(a | b) ≡ filter(b | a)  (as multisets)"""
    a, b = pred_a(), pred_b()
    left = _sort_rows(ann.filter(a | b).collect())
    right = _sort_rows(ann.filter(b | a).collect())
    assert left == right


# ---------------------------------------------------------------------------
# 7. Comparison reflexivity: filter(col == col).count() ≡ count()
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann())
def test_eq_reflexive(ann):
    """filter(col("x") == col("x")).count() == ann.count()"""
    assert ann.filter(col("x") == col("x")).count() == ann.count()


# ---------------------------------------------------------------------------
# 8. Comparison anti-reflexivity: filter(col != col).count() == 0
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann())
def test_ne_anti_reflexive(ann):
    """filter(col("x") != col("x")).count() == 0"""
    assert ann.filter(col("x") != col("x")).count() == 0


# ---------------------------------------------------------------------------
# 9. Select preserves row count: select(*cols).count() ≡ count()
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann())
def test_select_preserves_count(ann):
    """select(*cols).count() ≡ count()"""
    assert ann.select("x", "y").count() == ann.count()


# ---------------------------------------------------------------------------
# 10. head bound: head(n).count() <= min(n, count())
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann(), n=st.integers(min_value=0, max_value=30))
def test_head_bounded(ann, n):
    """head(n).count() <= min(n, ann.count())"""
    assert ann.head(n).count() <= min(n, ann.count())


# ---------------------------------------------------------------------------
# 11. Distinct idempotence: distinct().distinct().count() ≡ distinct().count()
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann())
def test_distinct_idempotent(ann):
    """distinct().distinct().count() ≡ distinct().count()"""
    once = ann.distinct().count()
    twice = ann.distinct().distinct().count()
    assert once == twice


# ---------------------------------------------------------------------------
# 12. Rename round-trip: rename(x->tmp).rename(tmp->x) ≡ ann (as multisets)
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann())
def test_rename_round_trip(ann):
    """rename(x->'_tmp').rename(_tmp->x) produces the same rows."""
    round_trip = _sort_rows(ann.rename(x="_tmp").rename(_tmp="x").collect())
    original = _sort_rows(ann.collect())
    assert round_trip == original


# ---------------------------------------------------------------------------
# 13. Arithmetic: with_columns(z=col(x)+col(y)) produces x+y for each row
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann())
def test_with_columns_arithmetic(ann):
    """with_columns(z=col("x")+col("y")) gives z == x+y for every row."""
    annotated = ann.with_columns(z=col("x") + col("y")).collect()
    for row in annotated:
        assert row["z"] == row["x"] + row["y"]


# ---------------------------------------------------------------------------
# 14. isin equivalence: col.isin([a,b]) ≡ (col==a)|(col==b)
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann())
def test_isin_equiv_to_disjunction(ann):
    """col("s").isin(["a","b"]) ≡ (col("s")=="a")|(col("s")=="b")"""
    via_isin = ann.filter(col("s").isin(["a", "b"])).count()
    via_or = ann.filter((col("s") == "a") | (col("s") == "b")).count()
    assert via_isin == via_or


# ---------------------------------------------------------------------------
# 15. gt implies ge: filter(x>0).count() <= filter(x>=0).count()
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann())
def test_gt_implies_ge(ann):
    """filter(x>0).count() <= filter(x>=0).count()"""
    gt = ann.filter(col("x") > 0).count()
    ge = ann.filter(col("x") >= 0).count()
    assert gt <= ge


# ---------------------------------------------------------------------------
# 16. collect/count consistency: count() == len(collect())
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann())
def test_count_matches_collect_len(ann):
    """ann.count() == len(ann.collect())"""
    assert ann.count() == len(ann.collect())


# ---------------------------------------------------------------------------
# 17. Filter De Morgan: count(filter(a))+count(filter(~a)) == count()
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann(), pred=_pred_strat)
def test_filter_de_morgan_partition(ann, pred):
    """filter(a).count() + filter(~a).count() == ann.count()"""
    p = pred()
    assert ann.filter(p).count() + ann.filter(~p).count() == ann.count()


# ---------------------------------------------------------------------------
# 18. And-identity with tautology: filter(a & (x==x)) ≡ filter(a)
# ---------------------------------------------------------------------------

@settings(**_SETTINGS)
@given(ann=_small_ann(), pred=_pred_strat)
def test_and_tautology_identity(ann, pred):
    """filter(a & (col("x")==col("x"))) ≡ filter(a) because (x==x) is always True"""
    p = pred()
    tautology = col("x") == col("x")
    assert (
        _sort_rows(ann.filter(p & tautology).collect())
        == _sort_rows(ann.filter(p).collect())
    )
