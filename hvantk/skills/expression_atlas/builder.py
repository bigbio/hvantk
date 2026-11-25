"""AnnData builder for the EBI Expression Atlas resource.

Owns the ``build_expression_atlas`` builder. Turns an Expression
Atlas baseline bulk-RNA-seq expression TSV plus its SDRF metadata file into
an ``ExpressionMatrix``: one observation per sequencing run (the default), or,
with ``observations="assay_group"``, one per assay group from the gene-level
``<accession>-tpms.tsv`` and ``<accession>-configuration.xml``. A matrix with
no non-zero finite value, or with an observation whose values are all zero, is
refused.

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
    quartiles) rather than per run, and then it cannot be built per run (it can per
    assay group, with ``observations="assay_group"``; see
    ``_resolve_assay_group_inputs``). So only the
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
    sdrf_path, accession = _find_sdrf(root)

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
            f"{accession}-transcripts-tpms.tsv) to --raw-dir, or build one "
            "observation per assay group with --plugin-arg observations=assay_group."
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


def _find_sdrf(root):
    """``(sdrf_path, accession)`` for the one ``*.condensed-sdrf.tsv`` under *root*."""
    sdrf_candidates = sorted(glob.glob(os.path.join(root, "*.condensed-sdrf.tsv")))
    if len(sdrf_candidates) != 1:
        raise FileNotFoundError(
            f"Expression Atlas builder: expected exactly one *.condensed-sdrf.tsv "
            f"under {root!r}, found {sdrf_candidates}. Point --raw-dir at a single "
            "accession's download directory."
        )
    sdrf_path = sdrf_candidates[0]
    return sdrf_path, _accession_of(sdrf_path)


def _accession_of(sdrf_path):
    """The accession an SDRF's file name (``<accession>.condensed-sdrf.tsv``) names."""
    name = os.path.basename(str(sdrf_path))
    suffix = ".condensed-sdrf.tsv"
    return name[: -len(suffix)] if name.endswith(suffix) else "<accession>"


def _resolve_assay_group_inputs(parsed_input):
    """``(expression_matrix, configuration, sdrf)`` paths for ``observations="assay_group"``.

    A mapping must carry ``expression_matrix``, ``configuration`` and ``sdrf`` keys. A
    raw dir must hold exactly one ``<accession>.condensed-sdrf.tsv`` and, for that
    accession, the gene-level ``<accession>-tpms.tsv`` (one column per assay group)
    and ``<accession>-configuration.xml`` (which runs make up each group); both come
    with an unfiltered ``hvantk download expression-atlas --accession <accession>``.
    """
    if isinstance(parsed_input, dict):
        missing = [
            k
            for k in ("expression_matrix", "configuration", "sdrf")
            if k not in parsed_input
        ]
        if missing:
            raise KeyError(
                "observations='assay_group' needs parsed_input keys expression_matrix, "
                f"configuration and sdrf; missing {missing}"
            )
        return (
            str(parsed_input["expression_matrix"]),
            str(parsed_input["configuration"]),
            str(parsed_input["sdrf"]),
        )

    root = str(parsed_input)
    sdrf_path, accession = _find_sdrf(root)
    matrix = os.path.join(root, f"{accession}-tpms.tsv")
    configuration = os.path.join(root, f"{accession}-configuration.xml")
    missing = [p for p in (matrix, configuration) if not os.path.isfile(p)]
    if missing:
        raise FileNotFoundError(
            f"Expression Atlas builder (observations=assay_group): {missing} not found "
            f"under {root!r}. The route reads the gene-level {accession}-tpms.tsv and "
            f"{accession}-configuration.xml; download them with `hvantk download "
            f"expression-atlas --accession {accession}` (no --config-path filter)."
        )
    logger.info(
        "Expression Atlas builder: one observation per assay group, from %s "
        "(groups from %s)",
        matrix,
        configuration,
    )
    return matrix, configuration, sdrf_path


def _zero_profile(X, block_rows=64):
    """``(any_nonzero, zero_rows)`` of dense *X*, scanned in row blocks.

    ``any_nonzero``: whether any value is finite and non-zero. ``zero_rows``: the rows
    that have finite values and every one of them 0. A row with no finite value at all
    (an all-missing sample) is neither.
    """
    import numpy as np

    any_nonzero = False
    zero_rows = []
    for start in range(0, X.shape[0], block_rows):
        block = np.asarray(X[start : start + block_rows])
        finite = np.isfinite(block)
        has_nonzero = (finite & (block != 0)).any(axis=1)
        any_nonzero = any_nonzero or bool(has_nonzero.any())
        zero_rows.extend(start + np.flatnonzero(finite.any(axis=1) & ~has_nonzero))
    return any_nonzero, zero_rows


def _require_expression(adata, expression_matrix_path, sdrf_path, observations):
    """Raise when the matrix, or any observation in it, is all zero.

    EBI's E-MTAB-6814-transcripts-tpms.tsv is all zero upstream (168,009 transcripts x
    287 runs). It is the accession's only per-run TPM file, so the resolver picks it,
    and the build wrote an all-zero ``.h5ad`` that read downstream as "not expressed".
    One all-zero run or group is the same defect on a smaller scale: a real TPM profile
    is never all zero. An observation with no finite value at all (an all-missing
    sample) stays allowed; that is "no data", not "not expressed".
    """
    any_nonzero, zero_rows = _zero_profile(adata.X)
    if any_nonzero and not zero_rows:
        return
    hint = ""
    if observations == "run":
        accession = _accession_of(sdrf_path)
        hint = (
            " Some Expression Atlas per-run files are all zero upstream (E-MTAB-6814's "
            "transcripts-tpms.tsv, for one). If the raw dir also holds the gene-level "
            f"{accession}-tpms.tsv (one column per assay group) and "
            f"{accession}-configuration.xml, build one observation per assay group from "
            "them instead: --plugin-arg observations=assay_group."
        )
    if not any_nonzero:
        raise ValueError(
            f"Expression Atlas builder: {expression_matrix_path} has no non-zero finite "
            f"value ({adata.n_obs} x {adata.n_vars} after parsing); refusing to write "
            f"an all-zero dataset.{hint}"
        )
    names = [str(adata.obs_names[i]) for i in zero_rows[:5]]
    raise ValueError(
        f"Expression Atlas builder: in {expression_matrix_path}, {len(zero_rows)} of "
        f"{adata.n_obs} observations have only zero values (e.g. {names}). A real "
        "expression profile is never all zero, so the file is broken for them; "
        f"refusing to write a dataset that reads them as 'not expressed'.{hint}"
    )


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
    observations: str = "run",
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
        already gene-level (see ``_resolve_expression_atlas_inputs``), and with
        ``observations="assay_group"``, whose input is always gene-level.
    observations : str
        ``"run"`` (default): one observation per sequencing run, from the TPM file
        whose columns are the SDRF's runs (``_resolve_expression_atlas_inputs``).
        ``"assay_group"``: one observation per assay group, from the gene-level
        ``<accession>-tpms.tsv`` and ``<accession>-configuration.xml``
        (``_resolve_assay_group_inputs``); ``X`` is the group's median TPM over its
        runs, and ``obs`` carries the group's label, runs and the SDRF values its
        runs share (``create_anndata_from_assay_groups``). An explicit opt-in, so a
        per-run build never silently becomes a per-group one; use it when the
        accession has no usable per-run file (E-MTAB-6814's is all zero upstream).
        Set with ``--plugin-arg observations=assay_group``. ``gene_column``,
        ``gene_name_column``, ``transcript_id_column`` and
        ``extra_annotation_columns`` apply to the per-run route only.
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
        ``var`` is keyed by the gene id like a gene-level build. With
        ``observations="assay_group"``, ``obs`` has one row per assay group
        (index name ``"assay_group"``) and ``var`` is keyed by the gene id.

    Raises
    ------
    ValueError
        If ``observations`` is not ``"run"`` or ``"assay_group"``, or the built
        matrix, or any observation in it, has finite values that are all zero
        (``_require_expression``): an all-zero upstream file must not become a
        dataset that reads as "not expressed".
    """
    from hvantk.core.models import ExpressionMatrix
    from hvantk.core.models.anndata_utils import annotate_column_summary_ad
    from hvantk.skills.expression_atlas.shared.expression_atlas import (
        convert_sdrf_to_dataframe,
        create_anndata_from_assay_groups,
        create_anndata_from_expression_atlas,
    )

    if observations not in ("run", "assay_group"):
        raise ValueError(
            f"observations={observations!r} is not supported; use 'run' (default, one "
            "observation per sequencing run) or 'assay_group' (one per assay group, "
            "the group's median TPM)."
        )

    if observations == "assay_group":
        if level not in ("transcript", "gene"):
            raise ValueError(
                f"level={level!r} is not supported; use 'transcript' or 'gene' (a "
                "no-op with observations='assay_group', whose input is gene-level)."
            )
        expression_matrix_path, configuration_file, sdrf_file = (
            _resolve_assay_group_inputs(parsed_input)
        )
        adata = create_anndata_from_assay_groups(
            expression_matrix_path,
            configuration_file,
            metadata_df=convert_sdrf_to_dataframe(sdrf_file),
            delimiter=delimiter,
        )
    else:
        expression_matrix_path, sdrf_file = _resolve_expression_atlas_inputs(
            parsed_input
        )
        adata = create_anndata_from_expression_atlas(
            expression_matrix_path=expression_matrix_path,
            metadata_df=convert_sdrf_to_dataframe(sdrf_file),
            gene_id_column=gene_column,
            gene_name_column=gene_name_column,
            transcript_id_column=transcript_id_column,
            level=level,
            delimiter=delimiter,
            extra_annotation_columns=extra_annotation_columns,
        )

    _require_expression(adata, expression_matrix_path, sdrf_file, observations)
    annotate_column_summary_ad(adata)

    return ExpressionMatrix.from_anndata(
        adata, provenance=ctx.provenance(schema_id="expression-atlas-dataset-v1")
    )
