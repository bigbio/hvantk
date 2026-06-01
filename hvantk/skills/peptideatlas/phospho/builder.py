"""Pandas builder for the PeptideAtlas Human Phospho resource.

PeptideAtlas does not have a dedicated Hail Table or AnnData representation
in hvantk today - the downstream consumer is
:func:`hvantk.algorithms.ptm.pipeline.ptm_build_pipeline`, which reads the intermediate
wide TSV produced by the dataset class directly (see
``peptideatlas_tsv`` in ``PTMBuildConfig``). This builder therefore delegates
to the existing parse code in
:mod:`hvantk.skills.peptideatlas.phospho.shared.datasets` and returns the
result as a :class:`pandas.DataFrame` so plugin consumers can validate the
schema and row contents the same way they do for any other tabular dataset.

The compound dataset key for the loader is ``peptideatlas:phospho``.
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


def build_peptideatlas_phospho(
    parsed_input,
    ctx,
    *,
    overwrite: bool = False,
):
    """Phase B builder — returns an AnnotationTable.

    ``parsed_input`` is whatever the plugin's parse_fn produced. For peptideatlas
    this is the path to the intermediate wide TSV (peptideatlas-phospho-*.tsv)
    that ``parse_raw_dir`` writes.

    Parameters
    ----------
    parsed_input : str | Path
        Path to the intermediate TSV file produced by parse_raw_dir.
    ctx : hvantk.core.models.BuildContext
        Platform-provided context. The plugin supplies schema_id via
        ``ctx.provenance(schema_id=...)``.
    overwrite : bool, optional
        Unused under the Phase B contract (the platform handles output writing).
        Accepted for backward compatibility only.

    Returns
    -------
    hvantk.core.models.AnnotationTable
        The intermediate phospho-site table wrapped with Provenance.
    """
    import pandas as pd
    from hvantk.core.models import AnnotationTable

    tsv_path = str(parsed_input)
    logger.info("Loading intermediate TSV %s", tsv_path)
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    return AnnotationTable.from_pandas(
        df, provenance=ctx.provenance(schema_id="peptideatlas-phospho-v1")
    )
