"""AnnData builder for UCSC Cell Browser collections.

Owns the Phase B ``build_ucsc_cellbrowser`` builder. Turns a UCSC Cell
Browser expression TSV plus its metadata file into an ``anndata.AnnData``
object (cells x genes), returned as an ``ExpressionMatrix`` for the platform
to persist.

Shared AnnData helpers (``build_anndata_metadata``,
``annotate_column_summary_ad``, ``save_anndata``) intentionally stay in
``hvantk/core/models/anndata_utils.py`` and ``hvantk/core/io/anndata_io.py``
because they are reused by every anndata builder (Expression Atlas, CPTAC,
...).
"""

from __future__ import annotations

import logging
import os

logger = logging.getLogger(__name__)

# Auto-select the backed builder for UCSC inputs larger than this threshold.
BACKED_BUILDER_THRESHOLD_BYTES = 1 * 1024 * 1024 * 1024  # 1 GiB

__all__ = [
    "build_ucsc_cellbrowser",
    "BACKED_BUILDER_THRESHOLD_BYTES",
]

# Map compound dataset names → Phase B schema IDs.
_SCHEMA_IDS: dict[str, str] = {
    "ucsc-cellbrowser:default":   "ucsc-cellbrowser-default-v1",
    "ucsc-cellbrowser:adult-ctx": "ucsc-cellbrowser-adult-ctx-v1",
    "ucsc-cellbrowser:dev-ctx":   "ucsc-cellbrowser-dev-ctx-v1",
}


def build_ucsc_cellbrowser(
    parsed_input,
    ctx,
    *,
    gene_column: str = "gene",
    delimiter: str = "\t",
    split_gene_field: bool = True,
    chunk_size: int = 500,
    backed: bool | None = None,
    column_batch: int = 64,
    backed_output_path: str | None = None,
    **params,
):
    """Phase B builder — returns an ExpressionMatrix.

    ``parsed_input`` must contain keys ``expression_matrix`` and ``metadata``.
    The per-dataset ``schema_id`` is resolved from ``ctx.dataset`` via
    ``_SCHEMA_IDS`` so that all three datasets (default, adult-ctx, dev-ctx)
    share one builder function while each stamps the correct schema.

    Backed vs in-memory dispatch:

    - In-memory (default for small inputs): loads the expression matrix
      directly into an ``AnnData`` object, annotates summary stats, and
      returns it. The platform's ``artifact.save()`` persists to disk.
    - Backed (auto-selected for inputs > ``BACKED_BUILDER_THRESHOLD_BYTES``,
      or forced via ``backed=True``): writes to disk while streaming. Phase B
      callers that need backed mode must pass ``backed_output_path`` (the
      final on-disk path) because backed writing requires materializing
      directly to a file. The platform then re-saves to the requested
      output via ``artifact.save()`` (no-op when paths match).
    """
    import anndata as ad

    from hvantk.core.models import ExpressionMatrix
    from hvantk.core.models.anndata_utils import (
        annotate_column_summary_ad,
        build_anndata_metadata,
    )
    from hvantk.skills.ucsc_cellbrowser.shared.ucsc import (
        build_ucsc_atlas_backed,
        coerce_obs_for_h5ad,
        create_anndata_from_ucsc_matrix,
        load_ucsc_metadata,
    )

    expression_matrix_path = str(parsed_input["expression_matrix"])
    metadata_path = str(parsed_input["metadata"])

    logger.info("Loading UCSC metadata from %s", metadata_path)
    metadata_df = load_ucsc_metadata(metadata_path, sep=delimiter)

    # Auto-select backed mode by input size unless caller forces.
    if backed is None:
        size = os.path.getsize(expression_matrix_path)
        backed = size > BACKED_BUILDER_THRESHOLD_BYTES
        logger.info(
            "build_ucsc_cellbrowser: auto-selected backed=%s "
            "(input size %.2f GiB, threshold %.2f GiB)",
            backed,
            size / (1024**3),
            BACKED_BUILDER_THRESHOLD_BYTES / (1024**3),
        )

    sid = _SCHEMA_IDS.get(ctx.dataset, "ucsc-cellbrowser-unknown-v1")

    if backed:
        if backed_output_path is None:
            raise ValueError(
                "Backed mode requires --plugin-arg backed_output_path=<path> "
                "(backed writes go directly to disk). Use the same path as "
                "--output."
            )
        provenance = {
            "hvantk_metadata": build_anndata_metadata("UCSC", expression_matrix_path)
        }
        build_ucsc_atlas_backed(
            expression_matrix_path=expression_matrix_path,
            output_path=backed_output_path,
            metadata_df=metadata_df,
            gene_column=gene_column,
            delimiter=delimiter,
            split_gene_field=split_gene_field,
            column_batch=column_batch,
            overwrite=False,
            uns=provenance,
        )
        # annotate_column_summary_ad would need to scan the full X matrix and
        # is intentionally skipped for backed atlases (v1 trade-off).
        logger.info(
            "Backed atlas built at %s; returning a backed AnnData handle.",
            backed_output_path,
        )
        adata = ad.read_h5ad(backed_output_path, backed="r")
        return ExpressionMatrix.from_anndata(
            adata, provenance=ctx.provenance(schema_id=sid)
        )

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

    return ExpressionMatrix.from_anndata(
        adata, provenance=ctx.provenance(schema_id=sid)
    )
