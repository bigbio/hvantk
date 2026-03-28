"""
Gene-level summary of QTL cascade evidence.

Aggregates variant-level cascade classifications to gene level and
overlays constraint metrics, disease-gene labels, and colocalization
posteriors.
"""

import logging
from typing import Optional

import pandas as pd

from hvantk.qtlcascade.constants import DEFAULT_COLOC_H4_THRESHOLD

logger = logging.getLogger(__name__)


def build_cascade_gene_summary(
    cascade_ht_path: str,
    output_path: str,
    constraint_ht_path: Optional[str] = None,
    disease_genes_ht_path: Optional[str] = None,
    coloc_df: Optional[pd.DataFrame] = None,
    overwrite: bool = False,
):
    """Aggregate cascade results to gene level.

    Parameters
    ----------
    cascade_ht_path : str
        Path to the cascade Hail Table (from ``build_cascade``).
    output_path : str
        Path to write the gene-level summary Hail Table.
    constraint_ht_path : str, optional
        Path to a gnomAD constraint table keyed by ``gene_id``.
        Expected field: ``oe_lof_upper`` (LOEUF).
    disease_genes_ht_path : str, optional
        Path to a disease-gene table keyed by ``gene_id`` (Ensembl gene ID).
        Must match the cascade table's gene_id key.  ClinGen/GenCC tables
        are typically keyed by gene_symbol — transform them before passing.
    coloc_df : pd.DataFrame, optional
        Coloc results with columns ``gene_id``, ``H4``, ``tissue``.
    overwrite : bool
        Overwrite existing output.

    Returns
    -------
    hl.Table
        Gene-level summary keyed by ``gene_id`` with fields:
        ``n_eqtl_variants``, ``n_pqtl_variants``, ``n_concordant``,
        ``n_discordant``, ``best_eqtl_pvalue``, ``best_pqtl_pvalue``,
        ``has_complete_cascade``, and optional overlay fields.
    """
    import hail as hl

    logger.info("Loading cascade table: %s", cascade_ht_path)
    ht = hl.read_table(cascade_ht_path)

    # Determine available fields
    row_fields = set(ht.row)

    # Build aggregation expressions
    agg_exprs = dict(
        n_eqtl_variants=hl.agg.count_where(hl.is_defined(ht.eqtl_beta)),
        n_pqtl_variants=hl.agg.count_where(hl.is_defined(ht.pqtl_beta)),
        n_concordant=hl.agg.count_where(ht.cascade_class == "eqtl_mediated"),
        n_discordant=hl.agg.count_where(ht.cascade_class == "discordant"),
        best_eqtl_pvalue=hl.agg.min(ht.eqtl_pvalue),
        best_pqtl_pvalue=hl.agg.min(ht.pqtl_pvalue),
    )

    if "tissue" in row_fields:
        agg_exprs["tissues"] = hl.agg.collect_as_set(ht.tissue)

    if "gene_symbol" in row_fields:
        agg_exprs["gene_symbol"] = hl.agg.take(ht.gene_symbol, 1)[0]

    logger.info("Aggregating cascade to gene level")
    gene_ht = ht.group_by(ht.gene_id).aggregate(**agg_exprs)
    gene_ht = gene_ht.key_by("gene_id")

    # Derived field
    gene_ht = gene_ht.annotate(
        has_complete_cascade=gene_ht.n_concordant > 0,
    )

    # ------------------------------------------------------------------
    # Optional overlays
    # ------------------------------------------------------------------

    if constraint_ht_path:
        logger.info("Overlaying constraint metrics from %s", constraint_ht_path)
        constraint_ht = hl.read_table(constraint_ht_path)
        constraint_idx = constraint_ht.index(gene_ht.gene_id)
        gene_ht = gene_ht.annotate(
            oe_lof_upper=constraint_idx.oe_lof_upper,
        )

    if disease_genes_ht_path:
        logger.info("Overlaying disease-gene labels from %s", disease_genes_ht_path)
        disease_ht = hl.read_table(disease_genes_ht_path)
        gene_ht = gene_ht.annotate(
            is_disease_gene=hl.is_defined(disease_ht.index(gene_ht.gene_id)),
        )

    if coloc_df is not None and not coloc_df.empty:
        logger.info("Overlaying colocalization posteriors")
        coloc_summary = (
            coloc_df.groupby("gene_id")
            .agg(
                coloc_max_h4=("H4", "max"),
                coloc_n_tissues=(
                    "H4",
                    lambda x: int((x > DEFAULT_COLOC_H4_THRESHOLD).sum()),
                ),
            )
            .reset_index()
        )
        coloc_ht = hl.Table.from_pandas(coloc_summary).key_by("gene_id")
        coloc_idx = coloc_ht.index(gene_ht.gene_id)
        gene_ht = gene_ht.annotate(
            coloc_max_h4=coloc_idx.coloc_max_h4,
            coloc_n_tissues=coloc_idx.coloc_n_tissues,
        )

    logger.info("Checkpointing gene summary to %s", output_path)
    return gene_ht.checkpoint(output_path, overwrite=overwrite)
