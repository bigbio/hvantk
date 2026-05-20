"""AnnData builder for UCSC Cell Browser collections.

This module owns ``build_ucsc_ad``, the canonical builder that turns a UCSC
Cell Browser expression TSV plus its metadata file into an ``anndata.AnnData``
object (cells x genes). It was migrated out of
``hvantk/tables/matrix_builders.py`` so that everything UCSC-specific
(builder, streaming/backed helpers, downloader, dataset class, tests,
fixtures, SKILL) lives under the plugin folder at
:mod:`hvantk.skills.ucsc_cellbrowser`.

The shared AnnData helpers (``build_anndata_metadata``,
``annotate_column_summary_ad``, ``save_anndata``) intentionally stay in
``hvantk/core/anndata_utils.py`` because they are reused by every anndata
builder (Expression Atlas, CPTAC, ...).
"""

from __future__ import annotations

import anndata as ad
import logging
import os
from typing import Optional

logger = logging.getLogger(__name__)

# Auto-select the backed builder for UCSC inputs larger than this threshold.
BACKED_BUILDER_THRESHOLD_BYTES = 1 * 1024 * 1024 * 1024  # 1 GiB

__all__ = [
    "build_ucsc_ad",
    "build_ucsc_cellbrowser",
    "BACKED_BUILDER_THRESHOLD_BYTES",
]

# Map compound dataset names → Phase B schema IDs.
_SCHEMA_IDS: dict[str, str] = {
    "ucsc-cellbrowser:default":   "ucsc-cellbrowser-default-v1",
    "ucsc-cellbrowser:adult-ctx": "ucsc-cellbrowser-adult-ctx-v1",
    "ucsc-cellbrowser:dev-ctx":   "ucsc-cellbrowser-dev-ctx-v1",
}


def build_ucsc_ad(
    expression_matrix_path: str,
    metadata_path: str,
    output_path: Optional[str] = None,
    gene_column: str = "gene",
    delimiter: str = "\t",
    split_gene_field: bool = True,
    overwrite: bool = False,
    chunk_size: int = 500,
    backed: bool | None = None,
    column_batch: int = 64,
) -> "ad.AnnData":
    """Build an AnnData object from UCSC Cell Browser expression + metadata.

    Parameters
    ----------
    expression_matrix_path : str
        Path to expression TSV (genes x cells).
    metadata_path : str
        Path to metadata TSV.
    output_path : str, optional
        If provided, save the AnnData as ``.h5ad``.
    gene_column : str
        Name of the gene identifier column (default ``"gene"``).
    delimiter : str
        Column delimiter (default tab).
    split_gene_field : bool
        Split pipe-separated gene names, keeping the first element.
    overwrite : bool
        Allow overwriting *output_path* if it exists.
    chunk_size : int
        Gene rows per streaming chunk (default 500).
    backed : bool, optional
        Force the backed-write builder. When ``None`` (default), auto-select
        based on the expression matrix file size (``> BACKED_BUILDER_THRESHOLD_BYTES``
        triggers backed mode). When ``True``, ``output_path`` is required.
    column_batch : int
        Gene-column batch size used by the backed builder (default 64).
        Peak RAM scales with ``n_cells × column_batch``.

    Returns
    -------
    ad.AnnData
        Expression AnnData with metadata in ``obs`` and provenance in ``uns``.
        In backed mode, returns a read-backed AnnData handle (``backed='r'``).
    """
    import anndata as ad

    from hvantk.skills.ucsc_cellbrowser.shared.ucsc import (
        load_ucsc_metadata,
        create_anndata_from_ucsc_matrix,
        build_ucsc_atlas_backed,
        coerce_obs_for_h5ad,
    )
    from hvantk.core.models.anndata_utils import (
        build_anndata_metadata,
        annotate_column_summary_ad,
        save_anndata,
    )

    logger.info("Loading UCSC metadata from %s", metadata_path)
    metadata_df = load_ucsc_metadata(metadata_path, sep=delimiter)

    # Auto-select backed mode by input size unless caller forces.
    if backed is None:
        size = os.path.getsize(expression_matrix_path)
        backed = size > BACKED_BUILDER_THRESHOLD_BYTES
        logger.info(
            "build_ucsc_ad: auto-selected backed=%s (input size %.2f GiB, threshold %.2f GiB)",
            backed, size / (1024**3), BACKED_BUILDER_THRESHOLD_BYTES / (1024**3),
        )

    if backed:
        if output_path is None:
            raise ValueError("backed=True requires output_path (writes directly to disk).")
        provenance = {"hvantk_metadata": build_anndata_metadata("UCSC", expression_matrix_path)}
        build_ucsc_atlas_backed(
            expression_matrix_path=expression_matrix_path,
            output_path=output_path,
            metadata_df=metadata_df,
            gene_column=gene_column,
            delimiter=delimiter,
            split_gene_field=split_gene_field,
            column_batch=column_batch,
            overwrite=overwrite,
            uns=provenance,
        )
        # Return a backed-mode handle — shape + obs/var without materializing X.
        # annotate_column_summary_ad would need to scan the full X matrix and
        # is intentionally skipped for backed atlases (v1 trade-off).
        logger.info(
            "Backed atlas built at %s; skipping annotate_column_summary_ad "
            "(would materialize X). Returning a backed AnnData handle.",
            output_path,
        )
        return ad.read_h5ad(output_path, backed="r")

    logger.info("Creating AnnData from UCSC expression matrix (in-memory)")
    adata = create_anndata_from_ucsc_matrix(
        expression_matrix_path=expression_matrix_path,
        metadata_df=metadata_df,
        gene_column=gene_column,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
        chunk_size=chunk_size,
    )
    adata.uns["hvantk_metadata"] = build_anndata_metadata(
        "UCSC", expression_matrix_path
    )
    annotate_column_summary_ad(adata)
    # Coerce object-dtype obs columns before any write_h5ad — anndata's
    # vlen-string HDF5 writer chokes on NaN mixed with strings. Matches
    # the invariant enforced by build_ucsc_atlas_backed.
    adata.obs = coerce_obs_for_h5ad(adata.obs)
    if output_path:
        save_anndata(adata, output_path, overwrite=overwrite)
    return adata


def build_ucsc_cellbrowser(
    parsed_input,
    ctx,
    *,
    gene_column: str = "gene",
    delimiter: str = "\t",
    split_gene_field: bool = True,
    chunk_size: int = 500,
    backed=None,
    column_batch: int = 64,
    **params,
):
    """Phase B builder — returns an ExpressionMatrix.

    Same dispatch as build_ucsc_ad (in-memory vs backed) but without writing.
    The platform's run_builder_for_spec calls artifact.save() to persist.

    ``parsed_input`` must contain keys ``expression_matrix`` and ``metadata``.
    The per-dataset ``schema_id`` is resolved from ``ctx.dataset`` via
    ``_SCHEMA_IDS`` so that all three datasets (default, adult-ctx, dev-ctx)
    share one builder function while each stamps the correct schema.
    """
    from hvantk.core.models import ExpressionMatrix

    expression_matrix_path = str(parsed_input["expression_matrix"])
    metadata_path = str(parsed_input["metadata"])

    # Delegate to the legacy in-memory or backed builder.
    # output_path=None means "don't save" — the platform calls artifact.save().
    adata = build_ucsc_ad(
        expression_matrix_path=expression_matrix_path,
        metadata_path=metadata_path,
        output_path=None,
        gene_column=gene_column,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
        overwrite=False,
        chunk_size=chunk_size,
        backed=backed,
        column_batch=column_batch,
    )

    sid = _SCHEMA_IDS.get(ctx.dataset, "ucsc-cellbrowser-unknown-v1")
    return ExpressionMatrix.from_anndata(
        adata, provenance=ctx.provenance(schema_id=sid)
    )

