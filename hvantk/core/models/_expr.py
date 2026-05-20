"""Backend-agnostic expression algebra.

Algorithms build expressions with col() and operator overloads; the
expression tree is then compiled to a backend-native predicate by
hvantk.core.models._compile. The tree is frozen + introspectable so the
same expression can be logged, cached, or run on either backend without
changing the algorithm.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class Expr:
    def __and__(self, other: "Expr") -> "BinOp":
        return BinOp("and", self, _lit(other))

    def __or__(self, other: "Expr") -> "BinOp":
        return BinOp("or", self, _lit(other))

    def __invert__(self) -> "UnaryOp":
        return UnaryOp("not", self)

    def __eq__(self, other: Any) -> "BinOp":  # type: ignore[override]
        return BinOp("eq", self, _lit(other))

    def __ne__(self, other: Any) -> "BinOp":  # type: ignore[override]
        return BinOp("ne", self, _lit(other))

    def __gt__(self, other: Any) -> "BinOp":
        return BinOp("gt", self, _lit(other))

    def __ge__(self, other: Any) -> "BinOp":
        return BinOp("ge", self, _lit(other))

    def __lt__(self, other: Any) -> "BinOp":
        return BinOp("lt", self, _lit(other))

    def __le__(self, other: Any) -> "BinOp":
        return BinOp("le", self, _lit(other))

    def __add__(self, other: Any) -> "BinOp":
        return BinOp("add", self, _lit(other))

    def __sub__(self, other: Any) -> "BinOp":
        return BinOp("sub", self, _lit(other))

    def __mul__(self, other: Any) -> "BinOp":
        return BinOp("mul", self, _lit(other))

    def __truediv__(self, other: Any) -> "BinOp":
        return BinOp("div", self, _lit(other))

    def __pow__(self, other: Any) -> "BinOp":
        return BinOp("pow", self, _lit(other))

    def isin(self, values) -> "CallOp":
        return CallOp("isin", self, (Literal(tuple(values)),))

    def is_null(self) -> "CallOp":
        return CallOp("is_null", self, ())

    def is_not_null(self) -> "CallOp":
        return CallOp("is_not_null", self, ())

    # Disable hashing — these nodes are equality-overridden for AST construction,
    # not for set membership. Use repr() if you need a string key.
    __hash__ = None  # type: ignore[assignment]


@dataclass(frozen=True, eq=False)
class Col(Expr):
    name: str


@dataclass(frozen=True, eq=False)
class Literal(Expr):
    value: Any


@dataclass(frozen=True, eq=False)
class BinOp(Expr):
    op: str
    left: Expr
    right: Expr


@dataclass(frozen=True, eq=False)
class UnaryOp(Expr):
    op: str
    operand: Expr


@dataclass(frozen=True, eq=False)
class CallOp(Expr):
    name: str
    receiver: Expr
    args: tuple[Expr, ...] = ()


def col(name: str) -> Col:
    return Col(name)


def _lit(v: Any) -> Expr:
    return v if isinstance(v, Expr) else Literal(v)
