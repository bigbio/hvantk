"""
QTL cascade join and classification.

Outer-joins eQTL and pQTL Hail Tables on ``(locus, alleles, gene_id)``
and assigns each variant-gene pair to a mechanistic class:

    eqtl_mediated — concordant eQTL + pQTL (same effect direction)
    discordant    — both present, opposite directions
    eqtl_only     — eQTL without pQTL evidence
    pqtl_only     — pQTL without eQTL evidence

The triple-key join on gene_id prevents false cascades from LD: a variant
that is an eQTL for gene A and pQTL for gene B at the same locus is not
a cascade.

References
----------
- Fang et al. (2025) — mechanistic categories (Fig 4C)
"""

import logging
from typing import Optional

from hvantk.core.qtl_constants import (
    DEFAULT_EQTL_P_THRESHOLD,
    DEFAULT_PQTL_P_THRESHOLD,
)

from hvantk.core.backends import Backend, algorithm

logger = logging.getLogger(__name__)


@algorithm(
    backends=[Backend.HAIL],
    input_format="table",
    output_format="table",
    key_fields=["locus", "alleles", "gene_id"],
)
def build_cascade(
    eqtl_ht_path: str,
    pqtl_ht_path: str,
    output_path: str,
    eqtl_p_threshold: float = DEFAULT_EQTL_P_THRESHOLD,
    pqtl_p_threshold: float = DEFAULT_PQTL_P_THRESHOLD,
    tissue: Optional[str] = None,
    overwrite: bool = False,
):
    """Build the QTL cascade by outer-joining eQTL and pQTL tables.

    Both input tables must be keyed by ``(locus, alleles, gene_id)``.

    Parameters
    ----------
    eqtl_ht_path : str
        Path to the eQTL Hail Table.
    pqtl_ht_path : str
        Path to the pQTL Hail Table.
    output_path : str
        Path to write the cascade Hail Table.
    eqtl_p_threshold : float
        P-value threshold for eQTL significance.
    pqtl_p_threshold : float
        P-value threshold for pQTL significance.
    tissue : str, optional
        Filter both tables to this tissue before joining.
    overwrite : bool
        Overwrite existing output.

    Returns
    -------
    hl.Table
        Cascade table keyed by ``(locus, alleles, gene_id)`` with fields:
        ``eqtl_beta``, ``eqtl_se``, ``eqtl_pvalue``, ``pqtl_beta``,
        ``pqtl_se``, ``pqtl_pvalue``, ``cascade_class``,
        ``attenuation_ratio``, ``tissue``.
    """
    import hail as hl

    logger.info("Loading eQTL table: %s", eqtl_ht_path)
    eqtl_ht = hl.read_table(eqtl_ht_path)
    logger.info("Loading pQTL table: %s", pqtl_ht_path)
    pqtl_ht = hl.read_table(pqtl_ht_path)

    # Filter by tissue
    if tissue:
        if "tissue" in list(eqtl_ht.row):
            eqtl_ht = eqtl_ht.filter(eqtl_ht.tissue == tissue)
        if "tissue" in list(pqtl_ht.row):
            pqtl_ht = pqtl_ht.filter(pqtl_ht.tissue == tissue)

    # Apply p-value thresholds
    if eqtl_p_threshold > 0:
        eqtl_ht = eqtl_ht.filter(eqtl_ht.p_value <= eqtl_p_threshold)
    if pqtl_p_threshold > 0:
        pqtl_ht = pqtl_ht.filter(pqtl_ht.p_value <= pqtl_p_threshold)

    # Warn about multi-tissue data without a tissue filter
    eqtl_has_tissue = "tissue" in list(eqtl_ht.row)
    pqtl_has_tissue = "tissue" in list(pqtl_ht.row)
    if not tissue and (eqtl_has_tissue or pqtl_has_tissue):
        logger.warning(
            "No tissue filter provided but input table(s) contain a 'tissue' "
            "field. The outer join on (locus, alleles, gene_id) may cross-"
            "multiply rows from different tissues. Consider providing --tissue."
        )

    # Select and prefix fields to avoid name collisions
    eqtl_ht = eqtl_ht.select(
        eqtl_beta=eqtl_ht.beta,
        eqtl_se=eqtl_ht.se,
        eqtl_pvalue=eqtl_ht.p_value,
    )
    pqtl_ht = pqtl_ht.select(
        pqtl_beta=pqtl_ht.beta,
        pqtl_se=pqtl_ht.se,
        pqtl_pvalue=pqtl_ht.p_value,
    )

    # Outer join on (locus, alleles, gene_id)
    logger.info("Building cascade (outer join)")
    ht = eqtl_ht.join(pqtl_ht, how="outer")

    # Classify
    has_eqtl = hl.is_defined(ht.eqtl_beta)
    has_pqtl = hl.is_defined(ht.pqtl_beta)
    same_dir = has_eqtl & has_pqtl & (ht.eqtl_beta * ht.pqtl_beta > 0)

    ht = ht.annotate(
        cascade_class=hl.case()
        .when(same_dir, "eqtl_mediated")
        .when(has_eqtl & has_pqtl, "discordant")
        .when(has_eqtl, "eqtl_only")
        .when(has_pqtl, "pqtl_only")
        .or_missing(),
        attenuation_ratio=hl.if_else(
            same_dir,
            1.0 - hl.abs(ht.pqtl_beta) / hl.abs(ht.eqtl_beta),
            hl.missing(hl.tfloat64),
        ),
        tissue=hl.literal(tissue) if tissue else hl.missing(hl.tstr),
    )

    logger.info("Checkpointing cascade to %s", output_path)
    return ht.checkpoint(output_path, overwrite=overwrite)
