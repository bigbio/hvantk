"""One-call gene-burden FET: clean MT -> assembled CohortManifest-conformant gene table."""
from __future__ import annotations

from hvantk.algorithms.burden.aggregate import (
    assert_clean_mt,
    count_2x2,
    variant_reductions,
)
from hvantk.algorithms.burden.fet import run_gene_burden_fet


def run_from_mt(
    mt,
    *,
    gene_col,
    route_col,
    arm_col,
    key,
    carrier_mode="het",
    score_field=None,
    mtc=None,
):
    assert_clean_mt(
        mt, gene_col=gene_col, route_col=route_col, arm_col=arm_col, key=key
    )
    counts = count_2x2(
        mt,
        gene_field=gene_col,
        route_field=route_col,
        arm_field=arm_col,
        carrier_mode=carrier_mode,
    )
    reductions = variant_reductions(
        mt,
        gene_field=gene_col,
        route_field=route_col,
        arm_field=arm_col,
        carrier_mode=carrier_mode,
        score_field=score_field,
    )
    return run_gene_burden_fet(counts, reductions, mtc=mtc)
