"""AnnData on-disk I/O. Moved from core/models/anndata_utils.py in Phase Q
to honor the intra-core rule: core/models doesn't perform I/O.

Metadata enrichment helpers (build_anndata_metadata, annotate_column_summary_ad)
stay in core/models/anndata_utils.py because they're pure data transformations.
"""

import logging
import os

import anndata as ad

logger = logging.getLogger(__name__)


def save_anndata(
    adata: ad.AnnData,
    path: str,
    overwrite: bool = True,
) -> None:
    """Save AnnData to ``.h5ad`` file.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data object to save.
    path : str
        Output file path (should end in ``.h5ad``).
    overwrite : bool, optional
        If False and *path* already exists, raise :class:`FileExistsError`.

    Raises
    ------
    FileExistsError
        If *overwrite* is False and *path* exists.
    """
    if not overwrite and os.path.exists(path):
        raise FileExistsError(f"File already exists: {path}")

    logger.info("Saving AnnData (%d obs x %d var) to %s", adata.n_obs, adata.n_vars, path)
    adata.write_h5ad(path)


def load_anndata(path: str) -> ad.AnnData:
    """Load AnnData from ``.h5ad`` file.

    Parameters
    ----------
    path : str
        Path to the ``.h5ad`` file.

    Returns
    -------
    ad.AnnData
        The loaded annotated data object.
    """
    logger.info("Loading AnnData from %s", path)
    return ad.read_h5ad(path)
