"""Adjusted-genotype (``adj``) annotation.

Ported from gnomad_methods (``gnomad.utils.annotations``), MIT licensed, at
https://github.com/broadinstitute/gnomad_methods -- see NOTICE below.

WHY PORTED RATHER THAN IMPORTED: hvantk used exactly one function from the `gnomad`
package, ``annotate_adj``, and it is ~15 lines of Hail expression with no gnomAD data
behind it -- the work is all Hail. Carrying the dependency for it cost 30 packages
(`hgvs`, `ga4gh-vrs`, `onnx`, `onnxruntime`, `skl2onnx`, `psycopg2`, `protobuf`,
`sympy`, ...), and, worse, gnomad 0.8.2 pins `jsonschema<4` transitively, which
conflicts with hvantk's own declared `jsonschema>=4.0` for plugin and feature-spec
validation. That conflict is what made `poetry.lock` unresolvable-in-place and left it
stale for ~28 commits.

The thresholds are gnomAD's published defaults and are kept verbatim, so `adj` here
means the same thing it means in a gnomAD callset.

NOTICE
------
Copyright (c) 2018 Broad Institute, licensed under the MIT License. This module is a
derivative work: the expression logic is unchanged, the LGT/LAD fallback is retained,
and only the surrounding packaging differs.
"""
from __future__ import annotations

from typing import Union

import hail as hl

__all__ = ["annotate_adj", "get_adj_expr"]

# gnomAD's published adj thresholds. Changing these changes what `adj` means, so they
# are module constants rather than scattered literals.
ADJ_GQ = 20
ADJ_DP = 10
ADJ_AB = 0.2
ADJ_HAPLOID_DP = 5


def get_adj_expr(
    gt_expr: hl.expr.CallExpression,
    gq_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
    dp_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
    ad_expr: hl.expr.ArrayNumericExpression,
    adj_gq: int = ADJ_GQ,
    adj_dp: int = ADJ_DP,
    adj_ab: float = ADJ_AB,
    haploid_adj_dp: int = ADJ_HAPLOID_DP,
) -> hl.expr.BooleanExpression:
    """Return the boolean ``adj`` filter expression for a genotype.

    A genotype is ``adj`` when all three hold:

    - ``GQ >= adj_gq``
    - ``DP >= adj_dp`` (or ``>= haploid_adj_dp`` when the call is haploid)
    - for heterozygous calls only, the allele balance of each called allele is
      ``>= adj_ab``. Het-ref is a special case: only the alt allele is checked, since
      the ref side of a het-ref call carries no information about the alt.

    Homozygous calls skip the allele-balance test entirely -- hence the
    ``when(~is_het(), True)`` short-circuit.
    """
    return (
        (gq_expr >= adj_gq)
        & hl.if_else(
            gt_expr.is_haploid(), dp_expr >= haploid_adj_dp, dp_expr >= adj_dp
        )
        & (
            hl.case()
            .when(~gt_expr.is_het(), True)
            .when(gt_expr.is_het_ref(), ad_expr[gt_expr[1]] / dp_expr >= adj_ab)
            .default(
                (ad_expr[gt_expr[0]] / dp_expr >= adj_ab)
                & (ad_expr[gt_expr[1]] / dp_expr >= adj_ab)
            )
        )
    )


def annotate_adj(
    mt: hl.MatrixTable,
    adj_gq: int = ADJ_GQ,
    adj_dp: int = ADJ_DP,
    adj_ab: float = ADJ_AB,
    haploid_adj_dp: int = ADJ_HAPLOID_DP,
) -> hl.MatrixTable:
    """Annotate entries with the boolean ``adj`` field (assumes diploid).

    Falls back to the local-allele fields ``LGT``/``LAD`` when the dense ``GT``/``AD``
    are absent, which is the shape a VDS carries before ``hl.vds.to_dense_mt``. hvantk's
    own caller densifies first and checks for GT/AD, so the fallback is normally inert --
    it is kept so this function behaves identically to the upstream one for any other
    caller.
    """
    gt_expr = mt.LGT if "GT" not in mt.entry and "LGT" in mt.entry else mt.GT
    ad_expr = mt.LAD if "AD" not in mt.entry and "LAD" in mt.entry else mt.AD

    return mt.annotate_entries(
        adj=get_adj_expr(
            gt_expr, mt.GQ, mt.DP, ad_expr, adj_gq, adj_dp, adj_ab, haploid_adj_dp
        )
    )
