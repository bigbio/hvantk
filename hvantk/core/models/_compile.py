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
    AggOp,
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


# ---------- hail compiler ----------

def _compile_to_hail_with_field_accessor(expr: Expr, field_accessor) -> "Any":
    """Internal: compile an Expr to Hail using a custom field accessor callable.

    ``field_accessor(name)`` must return a Hail expression for the named field.
    This indirection allows compiling against Tables (``ht["field"]``) as well
    as against MatrixTable col/row scopes (``mt.col["field"]`` /
    ``mt.row["field"]``).
    """
    import hail as hl

    def go(e: Expr):
        if isinstance(e, Col):
            return field_accessor(e.name)
        if isinstance(e, Literal):
            return hl.literal(e.value) if isinstance(e.value, (list, tuple, set)) else e.value
        if isinstance(e, UnaryOp) and e.op == "not":
            return ~go(e.operand)
        if isinstance(e, BinOp):
            left = go(e.left)
            right = go(e.right)
            return _HAIL_BINOPS[e.op](left, right)
        if isinstance(e, CallOp):
            receiver = go(e.receiver)
            if e.name == "isin":
                values = go(e.args[0])
                if not isinstance(values, list):
                    values = list(values)
                return hl.literal(values).contains(receiver)
            if e.name == "is_null":
                return hl.is_missing(receiver)
            if e.name == "is_not_null":
                return hl.is_defined(receiver)
            raise NotImplementedError(f"hail: unknown CallOp {e.name!r}")
        raise NotImplementedError(f"hail: unknown Expr node {type(e).__name__}")

    return go(expr)


def compile_to_hail(expr: Expr, ht: "Any") -> "Any":
    """Compile an Expr to a Hail expression resolved against `ht`.

    ``ht`` must be a Hail :class:`~hail.Table`; field references use
    ``ht[field_name]`` which is Table-native field access.  For MatrixTable
    col/row scopes, use :func:`compile_to_hail_mt_col` or
    :func:`compile_to_hail_mt_row` instead.
    """
    import hail as hl  # local import: avoid importing hail at module load

    def go(e: Expr):
        if isinstance(e, Col):
            return ht[e.name]
        if isinstance(e, Literal):
            return hl.literal(e.value) if isinstance(e.value, (list, tuple, set)) else e.value
        if isinstance(e, UnaryOp) and e.op == "not":
            return ~go(e.operand)
        if isinstance(e, BinOp):
            left = go(e.left)
            right = go(e.right)
            return _HAIL_BINOPS[e.op](left, right)
        if isinstance(e, CallOp):
            receiver = go(e.receiver)
            if e.name == "isin":
                values = go(e.args[0])
                # Hail's .contains() lives on ArrayExpression, not TupleExpression;
                # always pass a list to hl.literal so we get an ArrayExpression.
                if not isinstance(values, list):
                    values = list(values)
                return hl.literal(values).contains(receiver)
            if e.name == "is_null":
                return hl.is_missing(receiver)
            if e.name == "is_not_null":
                return hl.is_defined(receiver)
            raise NotImplementedError(f"hail: unknown CallOp {e.name!r}")
        raise NotImplementedError(f"hail: unknown Expr node {type(e).__name__}")

    return go(expr)


_HAIL_BINOPS = {
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


def compile_to_hail_mt_col(expr: Expr, mt: "Any") -> "Any":
    """Compile an Expr to a Hail expression scoped to a MatrixTable's columns.

    Field references ``col("name")`` are resolved as ``mt.col["name"]``, which
    returns an expression bound to the MatrixTable's col scope.  The result is
    suitable for use with ``mt.filter_cols()``.

    Parameters
    ----------
    expr:
        An Expr tree whose ``Col`` nodes name column (obs) fields.
    mt:
        A ``hail.MatrixTable`` instance.
    """
    return _compile_to_hail_with_field_accessor(expr, lambda name: mt.col[name])


def compile_to_hail_mt_row(expr: Expr, mt: "Any") -> "Any":
    """Compile an Expr to a Hail expression scoped to a MatrixTable's rows.

    Field references ``col("name")`` are resolved as ``mt.row["name"]``, which
    returns an expression bound to the MatrixTable's row scope.  The result is
    suitable for use with ``mt.filter_rows()``.

    Parameters
    ----------
    expr:
        An Expr tree whose ``Col`` nodes name row (var) fields.
    mt:
        A ``hail.MatrixTable`` instance.
    """
    return _compile_to_hail_with_field_accessor(expr, lambda name: mt.row[name])


# ---------- aggregation compilers ----------

def compile_agg_to_pandas(agg: AggOp, group_df: pd.DataFrame):
    if agg.name == "count":
        return len(group_df)
    series = compile_to_pandas(agg.argument, group_df)
    return getattr(series, agg.name)()


def compile_agg_to_hail(agg: AggOp, ht: "Any"):
    import hail as hl
    if agg.name == "count":
        return hl.agg.count()
    inner = compile_to_hail(agg.argument, ht)
    return getattr(hl.agg, agg.name)(inner)
