"""Tests for compile_to_pandas / compile_to_hail.

This is the load-bearing piece of Phase A: if the two compilers diverge,
every downstream algorithm silently breaks on backend swap.
"""
from __future__ import annotations

import pandas as pd
import pytest

from hvantk.core.models._compile import compile_to_pandas
from hvantk.core.models._expr import col


@pytest.fixture()
def df():
    return pd.DataFrame(
        {
            "gene_symbol": ["BRCA1", "BRCA2", "TP53"],
            "score": [0.7, 0.4, 0.9],
            "chrom": ["chr17", "chr13", "chr17"],
        }
    )


def test_pandas_eq_literal(df):
    mask = compile_to_pandas(col("gene_symbol") == "BRCA1", df)
    assert mask.tolist() == [True, False, False]


def test_pandas_gt_literal(df):
    mask = compile_to_pandas(col("score") > 0.5, df)
    assert mask.tolist() == [True, False, True]


def test_pandas_and_or(df):
    mask = compile_to_pandas(
        (col("score") > 0.5) & (col("chrom") == "chr17"), df
    )
    assert mask.tolist() == [True, False, True]

    mask = compile_to_pandas(
        (col("score") > 0.5) | (col("chrom") == "chr13"), df
    )
    assert mask.tolist() == [True, True, True]


def test_pandas_not(df):
    mask = compile_to_pandas(~(col("score") > 0.5), df)
    assert mask.tolist() == [False, True, False]


def test_pandas_arithmetic(df):
    series = compile_to_pandas(col("score") * 2, df)
    assert series.tolist() == [1.4, 0.8, 1.8]


def test_pandas_isin(df):
    mask = compile_to_pandas(col("gene_symbol").isin(["BRCA1", "TP53"]), df)
    assert mask.tolist() == [True, False, True]


def test_pandas_is_null():
    df = pd.DataFrame({"x": [1, None, 3]})
    mask = compile_to_pandas(col("x").is_null(), df)
    assert mask.tolist() == [False, True, False]


def test_pandas_unknown_column_raises(df):
    with pytest.raises(KeyError):
        compile_to_pandas(col("nonexistent") == 1, df)


# ---------------------------------------------------------------------------
# Hail compiler tests (Task 4)
# ---------------------------------------------------------------------------

from hvantk.core.models._compile import compile_to_hail


def _make_ht():
    """Build a small Hail Table mirroring the pandas fixture."""
    import hail as hl

    rows = [
        {"gene_symbol": "BRCA1", "score": 0.7, "chrom": "chr17"},
        {"gene_symbol": "BRCA2", "score": 0.4, "chrom": "chr13"},
        {"gene_symbol": "TP53",  "score": 0.9, "chrom": "chr17"},
    ]
    return hl.Table.parallelize(
        rows,
        hl.tstruct(gene_symbol=hl.tstr, score=hl.tfloat64, chrom=hl.tstr),
    )


@pytest.mark.hail
def test_hail_eq_literal():
    ht = _make_ht()
    expr = compile_to_hail(col("gene_symbol") == "BRCA1", ht)
    filtered = ht.filter(expr).collect()
    assert [r.gene_symbol for r in filtered] == ["BRCA1"]


@pytest.mark.hail
def test_hail_and_or():
    ht = _make_ht()
    expr = (col("score") > 0.5) & (col("chrom") == "chr17")
    filtered = ht.filter(compile_to_hail(expr, ht)).collect()
    assert {r.gene_symbol for r in filtered} == {"BRCA1", "TP53"}


@pytest.mark.hail
def test_hail_not():
    ht = _make_ht()
    expr = ~(col("score") > 0.5)
    filtered = ht.filter(compile_to_hail(expr, ht)).collect()
    assert [r.gene_symbol for r in filtered] == ["BRCA2"]


@pytest.mark.hail
def test_hail_isin():
    ht = _make_ht()
    expr = col("gene_symbol").isin(["BRCA1", "TP53"])
    filtered = ht.filter(compile_to_hail(expr, ht)).collect()
    assert {r.gene_symbol for r in filtered} == {"BRCA1", "TP53"}


@pytest.mark.hail
def test_hail_arithmetic():
    ht = _make_ht()
    annotated = ht.annotate(score_doubled=compile_to_hail(col("score") * 2, ht))
    values = sorted(r.score_doubled for r in annotated.collect())
    assert values == pytest.approx([0.8, 1.4, 1.8])


@pytest.mark.hail
@pytest.mark.parametrize(
    "expr_factory",
    [
        lambda: col("score") > 0.5,
        lambda: col("gene_symbol") == "BRCA1",
        lambda: (col("score") > 0.5) & (col("chrom") == "chr17"),
        lambda: ~(col("score") > 0.5),
        lambda: col("gene_symbol").isin(["BRCA1", "TP53"]),
    ],
)
def test_hail_pandas_parity(expr_factory):
    """The two compilers MUST agree on which rows survive a predicate."""
    ht = _make_ht()
    df = ht.to_pandas()

    hail_mask = ht.filter(compile_to_hail(expr_factory(), ht)).collect()
    hail_genes = sorted(r.gene_symbol for r in hail_mask)

    pandas_mask = compile_to_pandas(expr_factory(), df)
    pandas_genes = sorted(df[pandas_mask]["gene_symbol"].tolist())

    assert hail_genes == pandas_genes
