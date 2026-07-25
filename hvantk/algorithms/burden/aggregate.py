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


def _carrier_variant_expr(mt, carrier_mode: str):
    """Per-(variant,sample) 'is a carrier' expr for the per-variant reductions."""
    check_carrier_mode(carrier_mode)
    if carrier_mode == "hom":
        return mt.GT.is_hom_var()
    # het / chet / homs_chet all count a het genotype as a carrier at the variant level;
    # chet vs homs_chet only change the *gene-level* qualifies rule (qualifies_expr).
    return (
        mt.GT.is_het() | mt.GT.is_hom_var()
        if carrier_mode == "homs_chet"
        else mt.GT.is_het()
    )


def count_2x2(mt, *, gene_field, route_field, arm_field, carrier_mode):
    """Per-(gene,route) distinct-sample carrier counts by arm -> pandas DataFrame."""
    _require_hail()
    cmt = build_per_gene_carrier_mt(mt, gene_field=gene_field, route_field=route_field)
    # carry the arm annotation onto the carrier MT's columns (same sample keys).
    # mt.cols() must be captured once: two separate calls return two distinct Table
    # objects, and Hail rejects combining expressions from different source objects.
    mt_cols = mt.cols()
    arm = mt_cols.select(_is_case=mt_cols[arm_field])
    cmt = cmt.annotate_cols(_is_case=arm[cmt.col_key]._is_case)
    q = qualifies_expr(cmt, carrier_mode)
    cmt = cmt.annotate_rows(
        a=hl.agg.count_where(q & cmt._is_case),
        b=hl.agg.count_where(q & ~cmt._is_case),
    )
    arm_totals = mt.aggregate_cols(
        hl.struct(
            n_case=hl.agg.count_where(mt[arm_field]),
            n_control=hl.agg.count_where(~mt[arm_field]),
        )
    )
    pdf = cmt.rows().to_pandas()
    pdf = pdf.rename(columns={gene_field: "gene", route_field: "route"})
    pdf["n_case"] = arm_totals.n_case
    pdf["n_control"] = arm_totals.n_control
    return pdf[["gene", "route", "a", "b", "n_case", "n_control"]]


def variant_reductions(
    mt,
    *,
    gene_field,
    route_field,
    arm_field,
    carrier_mode,
    score_field: str | None = None,
):
    """Per-(gene,route) reduction inputs from per-variant arm carrier counts.

    Note: ``ctrl_freq`` (below) is the control-CARRIER frequency of a variant
    (distinct control carriers / control samples), not an allele frequency.
    """
    _require_hail()
    if _is_array(mt[route_field]):
        mt = mt.explode_rows(mt[route_field])
    carr = _carrier_variant_expr(mt, carrier_mode)
    n_ctrl = mt.aggregate_cols(hl.agg.count_where(~mt[arm_field]))
    score = (
        mt[score_field]
        if (score_field and score_field in mt.row)
        else hl.missing(hl.tfloat64)
    )
    mt = mt.annotate_rows(
        _n_case_carr=hl.agg.count_where(carr & mt[arm_field]),
        _n_ctrl_carr=hl.agg.count_where(carr & ~mt[arm_field]),
        _score=score,
    )
    # _ctrl_freq is the control-carrier frequency (carriers / control samples), not
    # an allele frequency -- see `variant_reductions` docstring.
    mt = mt.annotate_rows(_ctrl_freq=mt._n_ctrl_carr / n_ctrl if n_ctrl else 0.0)
    rows = mt.rows()
    grouped = rows.group_by(gene=rows[gene_field], route=rows[route_field]).aggregate(
        n_case_var=hl.agg.count_where(rows._n_case_carr > 0),
        conc_num=hl.agg.filter(rows._n_case_carr > 0, hl.agg.max(rows._n_case_carr)),
        conc_den=hl.agg.sum(rows._n_case_carr),
        n_case_private=hl.agg.count_where(
            (rows._n_case_carr > 0) & (rows._n_ctrl_carr == 0)
        ),
        score_sum=hl.agg.filter(rows._n_case_carr > 0, hl.agg.sum(rows._score)),
        score_n=hl.agg.count_where(
            (rows._n_case_carr > 0) & hl.is_defined(rows._score)
        ),
        drivers=hl.agg.filter(
            rows._n_case_carr > 0,
            hl.agg.collect(hl.struct(cc=rows._n_case_carr, ctrl_freq=rows._ctrl_freq)),
        ),
    )
    pdf = grouped.to_pandas()
    # to_pandas() renders array<struct> cells as plain Python lists of
    # hail.utils.struct.Struct. Struct already supports d["cc"] (its __getitem__
    # delegates to _get_field), so this conversion isn't strictly required for the
    # downstream pandas layer's exact access pattern -- it's kept anyway to make the
    # cross-layer contract (algorithms/burden -> fet.py) explicit plain dicts rather
    # than a Hail-specific type, decoupling that layer from Hail regardless of version.
    pdf["drivers"] = pdf["drivers"].apply(
        lambda cell: [{"cc": x["cc"], "ctrl_freq": x["ctrl_freq"]} for x in cell]
    )
    return pdf
