"""One-call gene-burden FET: clean MT -> assembled CohortManifest-conformant gene table."""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import hail as hl
    import pandas as pd

from hvantk.algorithms.burden.aggregate import (
    assert_clean_mt,
    count_2x2,
    variant_reductions,
)
from hvantk.algorithms.burden.fet import run_gene_burden_fet


def run_from_mt(
    mt: "hl.MatrixTable",
    *,
    gene_col: str,
    route_col: str,
    arm_col: str,
    key: str,
    carrier_mode: str = "het",
    score_field: str | None = None,
    mtc: str | None = None,
) -> "pd.DataFrame":
    """Run the gene-burden Fisher-exact op end-to-end on a clean MatrixTable.

    Pipeline: validate the clean-input contract (``assert_clean_mt``) -> per-
    (gene[,route]) 2x2 carrier counts and variant-level reduction inputs
    (``aggregate.py``, distributed Hail) -> Fisher-exact test + min-p route
    selection + architecture reductions (``fet.py``, local pandas). The
    result is one row per gene, taken at its min-p ("winning") route, and is
    CohortManifest-conformant; a ``p_adj`` column is added when *mtc* is set.

    Parameters
    ----------
    mt : hl.MatrixTable
        Clean, split bi-allelic cohort MatrixTable with genotypes (``GT``),
        a gene column, a route (variant-class) column, and a binary
        case/control column.
    gene_col : str
        Row field naming the gene (values live in the *key* key-space).
    route_col : str
        Row field naming the variant-class "route" (e.g. lof/missense);
        may be ``array<str>``, in which case a variant contributes to every
        route it is tagged with.
    arm_col : str
        Column (sample) field: boolean, True for cases, False for controls.
    key : str
        Gene key space *gene_col* values live in; one of
        ``checks.KEY_SPACES`` (``"gene_id"``, ``"hgnc_id"``, ``"symbol"``).
    carrier_mode : str
        Genotype-to-carrier rule; one of ``checks.CARRIER_MODES``
        (``"het"``, ``"hom"``, ``"chet"``, ``"homs_chet"``). Default ``"het"``.
    score_field : str, optional
        Row field with a per-variant prediction score (e.g. REVEL) reduced
        into ``mean_score_case``. ``None`` (default) disables the reduction.
    mtc : str, optional
        Multiple-testing-correction method (``"bonferroni"`` or ``"bh"``)
        applied to the gene-level min-p; ``None`` (default) skips it.

    Returns
    -------
    pd.DataFrame
        One row per gene (column ``"gene"``) with its winning ``route``,
        Fisher ``minp``/``odds_ratio``, and that route's architecture
        reductions (``n_case_var``, ``conc``, ``driver_af``,
        ``mean_score_case``, ``frac_case_private``), plus ``p_adj`` when
        *mtc* is set.
    """
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
