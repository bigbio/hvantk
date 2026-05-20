"""Compile Expr trees to backend-native operations.

Each compiler walks the tree once and emits the equivalent native
expression (pandas Series of bools for filters, scalar Series for value
expressions). The two compilers MUST stay semantically aligned —
property-based parity tests enforce this.
"""
from __future__ import annotations

from typing import Any

import pandas as pd

from hvantk.core.models._expr import (
    BinOp,
    CallOp,
    Col,
    Expr,
    Literal,
    UnaryOp,
)


# ---------- pandas compiler ----------

def compile_to_pandas(expr: Expr, df: pd.DataFrame) -> pd.Series:
    """Compile an Expr to a pandas Series (boolean for predicates, scalar for value exprs)."""
    if isinstance(expr, Col):
        if expr.name not in df.columns:
            raise KeyError(f"column not in DataFrame: {expr.name!r}")
        return df[expr.name]
    if isinstance(expr, Literal):
        return expr.value
    if isinstance(expr, UnaryOp) and expr.op == "not":
        operand = compile_to_pandas(expr.operand, df)
        return ~operand
    if isinstance(expr, BinOp):
        left = compile_to_pandas(expr.left, df)
        right = compile_to_pandas(expr.right, df)
        return _PANDAS_BINOPS[expr.op](left, right)
    if isinstance(expr, CallOp):
        receiver = compile_to_pandas(expr.receiver, df)
        if expr.name == "isin":
            values = compile_to_pandas(expr.args[0], df)
            return receiver.isin(values)
        if expr.name == "is_null":
            return receiver.isna()
        if expr.name == "is_not_null":
            return receiver.notna()
        raise NotImplementedError(f"pandas: unknown CallOp {expr.name!r}")
    raise NotImplementedError(f"pandas: unknown Expr node {type(expr).__name__}")


_PANDAS_BINOPS = {
    "and": lambda a, b: a & b,
    "or":  lambda a, b: a | b,
    "eq":  lambda a, b: a == b,
    "ne":  lambda a, b: a != b,
    "gt":  lambda a, b: a > b,
    "ge":  lambda a, b: a >= b,
    "lt":  lambda a, b: a < b,
    "le":  lambda a, b: a <= b,
    "add": lambda a, b: a + b,
    "sub": lambda a, b: a - b,
    "mul": lambda a, b: a * b,
    "div": lambda a, b: a / b,
    "pow": lambda a, b: a ** b,
}


# ---------- hail compiler (stub; populated in Task 4) ----------

def compile_to_hail(expr: Expr, ht: "Any") -> "Any":
    raise NotImplementedError("compile_to_hail lands in Task 4")
