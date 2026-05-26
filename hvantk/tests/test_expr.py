"""Tests for the Expr algebra: column refs, literals, binary/unary ops, composition."""
from __future__ import annotations

import pytest

from hvantk.core.models._expr import (
    BinOp,
    CallOp,
    Col,
    Literal,
    UnaryOp,
    col,
)


def test_col_returns_col_node():
    c = col("score")
    assert isinstance(c, Col)
    assert c.name == "score"


def test_eq_returns_binop_with_literal_right():
    expr = col("score") == 0.5
    assert isinstance(expr, BinOp)
    assert expr.op == "eq"
    assert isinstance(expr.left, Col)
    assert isinstance(expr.right, Literal)
    assert expr.right.value == 0.5


def test_and_or_compose():
    e = (col("score") > 0.5) & (col("chrom") == "chr1")
    assert isinstance(e, BinOp)
    assert e.op == "and"
    assert isinstance(e.left, BinOp) and e.left.op == "gt"
    assert isinstance(e.right, BinOp) and e.right.op == "eq"


def test_not_returns_unary_not():
    e = ~(col("score") > 0.5)
    assert isinstance(e, UnaryOp)
    assert e.op == "not"


def test_arithmetic_ops_compose():
    e = col("a") + col("b") * 2
    assert isinstance(e, BinOp)
    assert e.op == "add"
    assert isinstance(e.right, BinOp)
    assert e.right.op == "mul"


def test_isin_returns_callop():
    e = col("gene").isin(["BRCA1", "BRCA2"])
    assert isinstance(e, CallOp)
    assert e.name == "isin"
    assert isinstance(e.receiver, Col)
    assert isinstance(e.args[0], Literal)
    assert e.args[0].value == ("BRCA1", "BRCA2")


def test_is_null_no_args():
    e = col("score").is_null()
    assert isinstance(e, CallOp)
    assert e.name == "is_null"
    assert e.args == ()


def test_expr_nodes_are_frozen():
    c = col("x")
    with pytest.raises((AttributeError, Exception)):
        c.name = "y"
