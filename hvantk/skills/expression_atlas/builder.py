"""AnnData builder for the EBI Expression Atlas resource.

Owns the ``build_expression_atlas`` builder. Turns an Expression
Atlas baseline bulk-RNA-seq expression TSV plus its SDRF metadata file into
an ``ExpressionMatrix`` (samples x genes).

The shared ``annotate_column_summary_ad`` helper lives in
``hvantk/core/models/anndata_utils`` and is reused by every anndata
builder (UCSC, CPTAC, ...).
"""

from __future__ import annotations

import csv
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
    the downloader writes); exactly one must be present. Its filename encodes the
    accession (``<accession>.condensed-sdrf.tsv``), which scopes every candidate
    expression matrix: SDRF and TPM candidates used to be counted independently,
    so a ``--raw-dir`` holding one accession's SDRF and a *different* accession's
    TPM file (e.g. a reused directory, or two downloads merged by hand) silently
    passed and combined mismatched accessions. A TPM candidate not named for the
    SDRF's accession is ignored for counting purposes, and if none remain this
    raises naming the mismatch.

    **Among same-accession TPM files, the expression matrix is chosen by its
    columns, not its name.** A plain ``--plugin-arg accession=...`` download (no
    ``config_path`` filter) fetches every file in the experiment's FTP directory
    via ``ftp.nlst()`` (``hvantk/skills/expression_atlas/cli.py``), so the
    directory can hold both ``<accession>-tpms.tsv`` and
    ``<accession>-transcripts-tpms.tsv``. The gene-level file can have one column
    per assay group (E-MTAB-6798's: ``g1`` ... ``g94``, each cell a string of
    quartiles) rather than per run, and then it cannot be built. So only the
    header line of each candidate is read, and a file is kept when its columns
    include at least one run id from the SDRF. Exactly one kept file is used;
    of several, the gene-level one (``*-tpms.tsv``, NOT
    ``*-transcripts-tpms.tsv``), since a per-run gene file needs no aggregation.
    A transcript-level pick is transcript-keyed and needs ``--plugin-arg
    level=gene`` before the generic gene-level consumers (``hvantk expression
    markers``/``summarize``, the PTM constraint loader) can read it. The file
    actually chosen is logged at INFO. This raises when no candidate has a run
    column, or the choice is still ambiguous.
    """
    if isinstance(parsed_input, dict):
        return str(parsed_input["expression_matrix"]), str(parsed_input["sdrf"])

    from hvantk.skills.expression_atlas.shared.expression_atlas import (
        convert_sdrf_to_dataframe,
    )

    root = str(parsed_input)
    sdrf_candidates = sorted(glob.glob(os.path.join(root, "*.condensed-sdrf.tsv")))
    if len(sdrf_candidates) != 1:
        raise FileNotFoundError(
            f"Expression Atlas builder: expected exactly one *.condensed-sdrf.tsv "
            f"under {root!r}, found {sdrf_candidates}. Point --raw-dir at a single "
            "accession's download directory."
        )
    sdrf_path = sdrf_candidates[0]
    accession = os.path.basename(sdrf_path)[: -len(".condensed-sdrf.tsv")]

    all_tpm = sorted(glob.glob(os.path.join(root, "*-tpms.tsv")))
    accession_tpm = [
        f for f in all_tpm if os.path.basename(f).startswith(accession + "-")
    ]
    if not accession_tpm:
        raise FileNotFoundError(
            f"Expression Atlas builder: the SDRF {sdrf_path!r} is for accession "
            f"{accession!r}, but no *-tpms.tsv under {root!r} is named for that "
            f"accession (found {all_tpm} for other accessions, or none at all). "
            "The SDRF and expression matrix must be for the same accession -- "
            "point --raw-dir at a single accession's download directory."
        )

    run_ids = set(map(str, convert_sdrf_to_dataframe(sdrf_path).index))
    per_run = [f for f in accession_tpm if run_ids.intersection(_header_columns(f))]
    if not per_run:
        raise FileNotFoundError(
            f"Expression Atlas builder: none of the TPM files for accession "
            f"{accession!r} under {root!r} ({accession_tpm}) has a column named for "
            f"a run in {sdrf_path!r}. Their columns are not run ids (e.g. "
            "assay-group summaries g1...gN), so there is no per-run expression "
            "matrix to build from. Add the accession's per-run TPM file (e.g. "
            f"{accession}-transcripts-tpms.tsv) to --raw-dir."
        )

    gene_tpm = [f for f in per_run if not f.endswith("-transcripts-tpms.tsv")]
    if len(per_run) == 1:
        tpm_path = per_run[0]
    elif len(gene_tpm) == 1:
        tpm_path = gene_tpm[0]
    else:
        raise FileNotFoundError(
            f"Expression Atlas builder: could not find exactly one gene-level "
            f"*-tpms.tsv among the TPM files whose columns are the SDRF's runs "
            f"(found {per_run}) for accession {accession!r} under {root!r}. Point "
            "--raw-dir at a single accession's download directory, or call "
            "build_expression_atlas directly with parsed_input="
            "{'expression_matrix': ..., 'sdrf': ...}."
        )

    if tpm_path in gene_tpm:
        logger.info(
            "Expression Atlas builder: using gene-level expression matrix %s",
            tpm_path,
        )
    else:
        logger.info(
            "Expression Atlas builder: using transcript-level expression matrix "
            "%s (pass --plugin-arg level=gene for a gene-keyed result)",
            tpm_path,
        )
    return tpm_path, sdrf_path


def _header_columns(path):
    """The column names on the first line of the tab-separated file *path*."""
    with open(path, newline="", encoding="utf-8") as fh:
        return next(csv.reader(fh, delimiter="\t"), [])


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
        directory. The SDRF must be unambiguous; the expression matrix is scoped
        to the SDRF's accession and is the TPM file whose columns are the SDRF's
        runs (the gene-level one when both files qualify; see
        ``_resolve_expression_atlas_inputs``).
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    level : str
        ``"transcript"`` (default) or ``"gene"``. Forwarded verbatim to
        ``create_anndata_from_expression_atlas`` -- see there for what each
        value does. Only relevant when the resolved input is transcript-level
        (the only file with a column per run, as for E-MTAB-6798, whose
        gene-level file has a column per assay group): set
        ``--plugin-arg level=gene`` on `hvantk reprocess` to make that build
        usable by the generic gene-level consumers (e.g. `hvantk expression
        markers`/`summarize` or the PTM constraint loader, neither of which
        understands a transcript-keyed `var`). A no-op when the input is
        already gene-level (see ``_resolve_expression_atlas_inputs``).
    accession, config_path : str, optional
        Accepted and otherwise ignored by the build stage. ``hvantk reprocess``
        forwards every ``--plugin-arg`` to all three lifecycle stages (download,
        parse, build); this plugin's download stage requires ``accession`` or
        ``config_path`` to know which experiment to fetch (there is no single
        canonical "latest" Expression Atlas file), but the build stage needs
        neither -- ``_resolve_expression_atlas_inputs`` locates the two input
        files within ``--raw-dir``. Declared explicitly (not swallowed
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
