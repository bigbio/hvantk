"""AnnData builder for the EBI Expression Atlas resource.

Owns the ``build_expression_atlas`` builder. Turns an Expression
Atlas baseline bulk-RNA-seq expression TSV plus its SDRF metadata file into
an ``ExpressionMatrix`` (samples x genes).

The shared ``annotate_column_summary_ad`` helper lives in
``hvantk/core/models/anndata_utils`` and is reused by every anndata
builder (UCSC, CPTAC, ...).
"""

from __future__ import annotations

import glob
import logging
import os

logger = logging.getLogger(__name__)


def _resolve_expression_atlas_inputs(parsed_input):
    """Normalise the builder input into ``(expression_matrix_path, sdrf_path)``.

    Modelled on the UCSC Cell Browser plugin's ``_resolve_ucsc_inputs``. Accepts
    either:

    * a mapping with ``expression_matrix`` / ``sdrf`` keys (direct builder callers,
      e.g. the snapshot tests), or
    * a path to the raw download directory. This is the ``hvantk reprocess``
      contract: the dataset declares no ``lifecycle.parse`` stage, so reprocess
      passes the raw dir straight to the builder ("no parse stage declared: the
      builder consumes the raw dir directly", ``reprocess_cli``). Before this
      function existed, that crashed every ``reprocess expression-atlas:dataset``
      build with ``TypeError: string indices must be integers`` (issue #383).

    The SDRF is matched by ``*.condensed-sdrf.tsv`` (the only file of that shape
    the downloader writes) and the expression matrix by ``*-tpms.tsv``. The latter
    deliberately also matches ``*-transcripts-tpms.tsv`` -- a raw download
    directory can hold BOTH the gene-level and transcript-level TPM file for one
    accession, and there is no way to prefer one over the other from the filename
    alone, so finding more than one candidate raises rather than silently picking
    one (same reasoning as UCSC's "found in multiple locations" guard).
    """
    if isinstance(parsed_input, dict):
        return str(parsed_input["expression_matrix"]), str(parsed_input["sdrf"])

    root = str(parsed_input)
    sdrf_candidates = sorted(glob.glob(os.path.join(root, "*.condensed-sdrf.tsv")))
    tpm_candidates = sorted(glob.glob(os.path.join(root, "*-tpms.tsv")))

    if len(sdrf_candidates) != 1:
        raise FileNotFoundError(
            f"Expression Atlas builder: expected exactly one *.condensed-sdrf.tsv "
            f"under {root!r}, found {sdrf_candidates}. Point --raw-dir at a single "
            "accession's download directory."
        )
    if len(tpm_candidates) != 1:
        raise FileNotFoundError(
            f"Expression Atlas builder: expected exactly one *-tpms.tsv under "
            f"{root!r}, found {tpm_candidates}. A download directory can hold both "
            "a gene-level <accession>-tpms.tsv and a transcript-level "
            "<accession>-transcripts-tpms.tsv; move the one you do not want to "
            "build out of --raw-dir, or call build_expression_atlas directly with "
            "parsed_input={'expression_matrix': ..., 'sdrf': ...}."
        )
    return tpm_candidates[0], sdrf_candidates[0]


def build_expression_atlas(
    parsed_input,
    ctx,
    *,
    gene_column: str = "Gene ID",
    gene_name_column: str = "Gene Name",
    transcript_id_column: str = "GeneID",
    level: str = "transcript",
    delimiter: str = "\t",
    extra_annotation_columns: tuple[str, ...] | list[str] | None = None,
    accession: str | None = None,
    config_path: str | None = None,
):
    """Plugin builder — returns an ExpressionMatrix.

    Parameters
    ----------
    parsed_input : dict or str
        Either a mapping with ``expression_matrix`` (path to gene x sample TSV)
        and ``sdrf`` (path to SDRF metadata file) keys, or a raw download
        directory containing exactly one of each (see
        ``_resolve_expression_atlas_inputs``).
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    level : str
        ``"transcript"`` (default) or ``"gene"``. Forwarded verbatim to
        ``create_anndata_from_expression_atlas`` -- see there for what each
        value does. Set via ``--plugin-arg level=gene`` on `hvantk reprocess`
        when a gene-keyed matrix is needed (e.g. for `hvantk expression
        markers`/`summarize` or the PTM constraint loader, neither of which
        understands a transcript-keyed `var`).
    accession, config_path : str, optional
        Accepted and otherwise ignored by the build stage. ``hvantk reprocess``
        forwards every ``--plugin-arg`` to all three lifecycle stages (download,
        parse, build); this plugin's download stage requires ``accession`` or
        ``config_path`` to know which experiment to fetch (there is no single
        canonical "latest" Expression Atlas file), but the build stage needs
        neither -- ``_resolve_expression_atlas_inputs`` locates the two input
        files by name within ``--raw-dir``. Declared explicitly (not swallowed
        by a bare ``**params``) so an actual typo in a ``--plugin-arg`` name
        still raises ``TypeError`` instead of being silently absorbed.

    Returns
    -------
    hvantk.core.models.ExpressionMatrix
        AnnData-backed ExpressionMatrix wrapped with Provenance. For a
        transcript-level export with the default ``level="transcript"``,
        ``var`` is keyed by the transcript id (see
        ``create_anndata_from_expression_atlas``); ``gene_column`` is kept as a
        ``var`` column rather than being dropped. With ``level="gene"``, a
        transcript-level export is summed onto the gene id instead, and
        ``var`` is keyed by the gene id like a gene-level build.
    """
    from hvantk.core.models import ExpressionMatrix
    from hvantk.core.models.anndata_utils import annotate_column_summary_ad
    from hvantk.skills.expression_atlas.shared.expression_atlas import (
        convert_sdrf_to_dataframe,
        create_anndata_from_expression_atlas,
    )

    expression_matrix_path, sdrf_file = _resolve_expression_atlas_inputs(parsed_input)

    metadata_df = convert_sdrf_to_dataframe(sdrf_file)
    adata = create_anndata_from_expression_atlas(
        expression_matrix_path=expression_matrix_path,
        metadata_df=metadata_df,
        gene_id_column=gene_column,
        gene_name_column=gene_name_column,
        transcript_id_column=transcript_id_column,
        level=level,
        delimiter=delimiter,
        extra_annotation_columns=extra_annotation_columns,
    )

    annotate_column_summary_ad(adata)

    return ExpressionMatrix.from_anndata(
        adata, provenance=ctx.provenance(schema_id="expression-atlas-dataset-v1")
    )
