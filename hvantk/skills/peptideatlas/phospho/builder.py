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
import os
from typing import Optional

logger = logging.getLogger(__name__)


def build_peptideatlas_phospho_tb(
    input_path: str,
    output_path: Optional[str] = None,
    overwrite: bool = False,
    **kwargs,
):
    """Build the PeptideAtlas Human Phospho intermediate table.

    ``input_path`` may point at either:

    * an ``atlas_build_*.tsv.zip`` archive (the raw upstream artefact), in
      which case the function parses it via
      :func:`hvantk.skills.peptideatlas.phospho.shared.datasets.parse_peptideatlas_zip`
      and writes the wide intermediate TSV via
      :func:`hvantk.skills.peptideatlas.phospho.shared.datasets.write_intermediate_tsv`;
      OR
    * an already-parsed ``peptideatlas-phospho-*.tsv`` file, in which case
      the function loads it directly with pandas.

    Parameters
    ----------
    input_path : str
        Path to the upstream zip OR a pre-parsed intermediate TSV.
    output_path : str, optional
        If provided and ``input_path`` is a zip, write the intermediate TSV
        here. Ignored when ``input_path`` is already a TSV.
    overwrite : bool
        If True, overwrite ``output_path`` when it already exists.
    **kwargs
        Reserved for future builder keyword arguments.

    Returns
    -------
    pandas.DataFrame
        The intermediate phospho-site table with the columns documented in
        ``shared.datasets._TSV_COLUMNS`` plus ``source_db`` and
        ``evidence_type``.
    """
    import pandas as pd

    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        parse_peptideatlas_zip,
        write_intermediate_tsv,
    )

    if input_path.endswith(".zip"):
        if output_path is None:
            raise ValueError(
                "output_path is required when input_path points at a "
                "PeptideAtlas atlas_build_*.tsv.zip archive."
            )
        if os.path.exists(output_path) and not overwrite:
            logger.info(
                "Intermediate TSV already exists, loading without re-parsing: %s",
                output_path,
            )
        else:
            logger.info("Parsing PeptideAtlas zip %s", input_path)
            sites = parse_peptideatlas_zip(input_path)
            write_intermediate_tsv(sites, output_path)
        tsv_path = output_path
    else:
        tsv_path = input_path

    logger.info("Loading intermediate TSV %s", tsv_path)
    return pd.read_csv(tsv_path, sep="\t", dtype=str)


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
