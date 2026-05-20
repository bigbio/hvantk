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
