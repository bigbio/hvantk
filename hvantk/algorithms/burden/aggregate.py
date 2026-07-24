"""Hail aggregation for gene-burden: clean MT -> per-gene[,route] carrier structures.

``build_per_gene_carrier_mt`` is the SHARED primitive: it collapses variants to a
per-(gene[,route]) x per-sample carrier MatrixTable with ``hets``/``homs``/``multi_het``
entries. The gene-set enrichment path (enrichex) is rewired onto it in a later task, so
its output schema deliberately matches the function it replaces
(``hvantk.algorithms.enrichex.burden.compute_per_gene_burden_mt``) when ``route_field``
is None.

Everything statistical happens later in ``fet.py`` on small collected frames; this module
only does the distributed per-sample counting and the fail-loud input contract.
"""
from __future__ import annotations

import logging

try:
    import hail as hl
except ModuleNotFoundError as exc:  # pragma: no cover - depends on env
    hl = None  # type: ignore
    _HAIL_IMPORT_ERROR = exc
else:
    _HAIL_IMPORT_ERROR = None

from hvantk.algorithms.burden.checks import (
    check_required_fields,
    check_key_space,
    check_carrier_mode,
)

logger = logging.getLogger(__name__)


def _require_hail() -> None:
    if hl is None:  # pragma: no cover - depends on env
        raise ImportError(
            "Hail is required for gene-burden aggregation."
        ) from _HAIL_IMPORT_ERROR


def _is_array(expr) -> bool:
    return str(expr.dtype).startswith("array<")


def build_per_gene_carrier_mt(mt, *, gene_field: str, route_field: str | None = None):
    """variants -> per-(gene[,route]) x per-sample carrier MT (hets/homs/multi_het).

    Identical to enrichex's ``compute_per_gene_burden_mt`` when ``route_field`` is
    None (that function is rewired onto this primitive in a later task).
    """
    _require_hail()
    if route_field is not None and _is_array(mt[route_field]):
        mt = mt.explode_rows(mt[route_field])

    mt = mt.filter_rows(hl.is_defined(mt[gene_field]))
    if route_field is not None:
        mt = mt.filter_rows(hl.is_defined(mt[route_field]))
        grouped = mt.group_rows_by(mt[gene_field], mt[route_field])
    else:
        grouped = mt.group_rows_by(mt[gene_field])

    return grouped.aggregate(
        hets=hl.agg.count_where(mt.GT.is_het()),
        homs=hl.agg.count_where(mt.GT.is_hom_var()),
        multi_het=hl.agg.count_where(mt.GT.is_het()) >= 2,
    )


def qualifies_expr(carrier_mt, carrier_mode: str):
    """Boolean entry expr: sample carries >=1 qualifying variant in this gene[,route]."""
    check_carrier_mode(carrier_mode)
    if carrier_mode == "het":
        return carrier_mt.hets > 0
    if carrier_mode == "hom":
        return carrier_mt.homs > 0
    if carrier_mode == "chet":
        return carrier_mt.multi_het
    return carrier_mt.multi_het | (carrier_mt.homs > 0)  # homs_chet


def assert_clean_mt(
    mt, *, gene_col: str, route_col: str, arm_col: str, key: str
) -> None:
    """Fail loud if the MT violates the clean-input contract.

    Runs the pure field-name checks first (cheap, no Spark job), then the
    data-dependent checks: bi-allelic rows and a defined, binary arm column.
    """
    _require_hail()
    check_key_space(key)
    check_required_fields(
        row_fields=list(mt.row),
        entry_fields=list(mt.entry),
        col_fields=list(mt.col),
        gene_col=gene_col,
        route_col=route_col,
        arm_col=arm_col,
    )

    n_multi = mt.filter_rows(hl.len(mt.alleles) != 2).count_rows()
    if n_multi:
        raise ValueError(
            f"{n_multi} multiallelic row(s): the burden op requires split, bi-allelic "
            "variants. Run hl.split_multi upstream."
        )

    stats = mt.aggregate_cols(
        hl.struct(
            n=hl.agg.count(),
            n_defined=hl.agg.count_where(hl.is_defined(mt[arm_col])),
            n_true=hl.agg.count_where(mt[arm_col]),
        )
    )
    if stats.n_defined != stats.n:
        raise ValueError(
            f"case/control column {arm_col!r} must be a defined boolean for every sample; "
            f"{stats.n - stats.n_defined}/{stats.n} are missing (sample QC is upstream)."
        )
    if stats.n_true == 0 or stats.n_true == stats.n:
        raise ValueError(
            f"case/control column {arm_col!r} must be binary: found "
            f"{stats.n_true} cases / {stats.n - stats.n_true} controls."
        )
