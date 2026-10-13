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
    the downloader writes); exactly one must be present.

    **The expression matrix prefers the gene-level file.** A plain
    ``--plugin-arg accession=...`` download (no ``config_path`` filter) fetches
    every file in the experiment's FTP directory via ``ftp.nlst()``
    (``hvantk/skills/expression_atlas/cli.py``), and an experiment that ships a
    transcript-level TPM file usually ships the gene-level ``<accession>-tpms.tsv``
    alongside it -- so requiring exactly one ``*-tpms.tsv`` match (as this function
    first did) could never build from an unfiltered download. When exactly one
    gene-level file (``*-tpms.tsv``, NOT ``*-transcripts-tpms.tsv``) is present, it
    is used and any transcript-level file is ignored. Otherwise, when exactly one
    ``*-transcripts-tpms.tsv`` is present (the experiment shipped only a
    transcript-level matrix, or a ``config_path`` download filtered to just that
    file), that one is used instead -- the result is transcript-keyed and needs
    ``--plugin-arg level=gene`` before the generic gene-level consumers (``hvantk
    expression markers``/``summarize``, the PTM constraint loader) can read it. The
    file actually chosen is logged at INFO. This raises only when the choice is
    still ambiguous (e.g. more than one accession's files mixed into one
    directory) or nothing matches.
    """
    if isinstance(parsed_input, dict):
        return str(parsed_input["expression_matrix"]), str(parsed_input["sdrf"])

    root = str(parsed_input)
    sdrf_candidates = sorted(glob.glob(os.path.join(root, "*.condensed-sdrf.tsv")))
    if len(sdrf_candidates) != 1:
        raise FileNotFoundError(
            f"Expression Atlas builder: expected exactly one *.condensed-sdrf.tsv "
            f"under {root!r}, found {sdrf_candidates}. Point --raw-dir at a single "
            "accession's download directory."
        )

    all_tpm = sorted(glob.glob(os.path.join(root, "*-tpms.tsv")))
    transcript_tpm = [f for f in all_tpm if f.endswith("-transcripts-tpms.tsv")]
    gene_tpm = [f for f in all_tpm if f not in transcript_tpm]

    if len(gene_tpm) == 1:
        tpm_path = gene_tpm[0]
        logger.info(
            "Expression Atlas builder: using gene-level expression matrix %s",
            tpm_path,
        )
    elif len(transcript_tpm) == 1:
        tpm_path = transcript_tpm[0]
        logger.info(
            "Expression Atlas builder: using transcript-level expression matrix "
            "%s (pass --plugin-arg level=gene for a gene-keyed result)",
            tpm_path,
        )
    else:
        raise FileNotFoundError(
            f"Expression Atlas builder: could not find exactly one gene-level "
            f"*-tpms.tsv or exactly one transcript-level *-transcripts-tpms.tsv "
            f"under {root!r}; found gene-level candidates {gene_tpm} and "
            f"transcript-level candidates {transcript_tpm}. Point --raw-dir at a "
            "single accession's download directory, or call build_expression_atlas "
            "directly with parsed_input={'expression_matrix': ..., 'sdrf': ...}."
        )
    return tpm_path, sdrf_candidates[0]


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
        directory. The SDRF must be unambiguous; the expression matrix prefers
        the experiment's gene-level TPM file when exactly one is present, falling
        back to the transcript-level one otherwise (see
        ``_resolve_expression_atlas_inputs``).
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    level : str
        ``"transcript"`` (default) or ``"gene"``. Forwarded verbatim to
        ``create_anndata_from_expression_atlas`` -- see there for what each
        value does. Only relevant when the resolved input is transcript-level
        (the experiment shipped no gene-level file, or a ``config_path``
        download filtered to just the transcript-level one): set
        ``--plugin-arg level=gene`` on `hvantk reprocess` to make that build
        usable by the generic gene-level consumers (e.g. `hvantk expression
        markers`/`summarize` or the PTM constraint loader, neither of which
        understands a transcript-keyed `var`). A no-op when the input is
        already gene-level, which is the common case for an unfiltered
        download (see ``_resolve_expression_atlas_inputs``).
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
