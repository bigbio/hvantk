"""
CPTAC dataset handling.

This module provides functions and dataset classes for downloading and
parsing data from CPTAC (Clinical Proteomic Tumor Analysis Consortium) via
the ``cptac`` Python package: phosphoproteomics (:class:`CPTACPhosphoDataset`)
and whole-proteome protein expression (:class:`CPTACExpressionDataset`). Both
share the cancer-type validation, clinical/tissue-type tagging, skip-and-continue
multi-cancer download loop, and raw-dir/cancer_type builder-input resolution
defined below.

Example usage:
    dataset = CPTACPhosphoDataset(cancer_type="brca")
    dataset.download("/data/cptac")

    sites = extract_phospho_sites(phospho_df, cancer_type="brca")
    write_intermediate_tsv(sites, "/data/cptac/cptac_phospho.tsv")
"""

import csv
import logging
import os
import re
import shutil
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from hvantk.core.plugin.drift_runner import _current_umask
from hvantk.skills.cptac.shared.constants import (
    CPTAC_CANCER_TYPES,
    CPTAC_CANCER_CLASS_MAP,
)

logger = logging.getLogger(__name__)

# Phospho amino acid descriptions
_PHOSPHO_AA_DESC = {
    "S": "Phosphoserine",
    "T": "Phosphothreonine",
    "Y": "Phosphotyrosine",
}

# Regex for parsing phospho site tokens like S65, T185, Y243
_SITE_RE = re.compile(r"([STY])(\d+)")

# Output TSV column order for intermediate file
_TSV_COLUMNS = [
    "accession",
    "gene_symbol",
    "position",
    "description",
    "amino_acid",
    "ensembl_xrefs",
    "sequence_length",
    "n_observations",
    "source_db",
    "evidence_type",
    "tissue_type",
    "cancer_type",
    "mean_intensity",
]


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------


def parse_phospho_site(site_str: str) -> List[Tuple[str, int]]:
    """Parse a CPTAC site string into (amino_acid, position) tuples.

    Accepts single sites (``"S65"``) and underscore-delimited multi-sites
    (``"T185_Y187"``).  Only S, T, and Y residues are recognised.

    Returns an empty list for empty strings or unrecognised residues.
    """
    if not site_str:
        return []
    return [(m.group(1), int(m.group(2))) for m in _SITE_RE.finditer(site_str)]


def extract_phospho_sites(
    phospho_df,
    cancer_type: str,
    tissue_type: str = "tumor",
) -> List[dict]:
    """Extract phospho sites from a CPTAC multi-index DataFrame.

    Parameters
    ----------
    phospho_df : pd.DataFrame
        DataFrame with a :class:`~pandas.MultiIndex` on columns whose first
        two levels are *Gene* and *Site*.
    cancer_type : str
        CPTAC cancer type identifier (e.g. ``"brca"``).
    tissue_type : str
        Tissue type label (``"tumor"`` or ``"normal"``).

    Returns
    -------
    list[dict]
        One dict per unique (gene, amino_acid, position) combination.
    """
    # Deduplicate: aggregate across columns that map to the same site
    agg: Dict[tuple, dict] = {}

    for col in phospho_df.columns:
        gene = col[0]
        site_str = col[1]
        db_id = col[3] if len(col) > 3 else ""

        parsed = parse_phospho_site(site_str)
        if not parsed:
            continue

        values = phospho_df[col].dropna()
        n_obs = len(values)
        mean_val = float(values.mean()) if n_obs > 0 else 0.0

        for aa, pos in parsed:
            key = (gene, aa, pos)
            if key in agg:
                existing = agg[key]
                total_obs = existing["n_observations"] + n_obs
                if total_obs > 0:
                    existing["mean_intensity"] = (
                        existing["mean_intensity"] * existing["n_observations"]
                        + mean_val * n_obs
                    ) / total_obs
                existing["n_observations"] = total_obs
            else:
                agg[key] = {
                    "accession": db_id,
                    "gene_symbol": gene,
                    "position": pos,
                    "description": _PHOSPHO_AA_DESC.get(aa, ""),
                    "amino_acid": aa,
                    "ensembl_xrefs": "",
                    "sequence_length": "",
                    "n_observations": n_obs,
                    "mean_intensity": round(mean_val, 4),
                    "source_db": "CPTAC",
                    "evidence_type": "mass_spectrometry",
                    "tissue_type": tissue_type,
                    "cancer_type": cancer_type,
                }

    return list(agg.values())


def write_intermediate_tsv(sites: List[dict], output_path: str) -> None:
    """Write extracted phospho sites to a 10-column TSV file."""
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=_TSV_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(sites)
    logger.info("Wrote %d phospho sites to %s", len(sites), output_path)


def write_matrix_csv(phospho_df, output_path: str) -> None:
    """Write phospho DataFrame as a matrix CSV (sites x samples).

    Multi-index columns are flattened to ``Gene_AminoAcid Position`` labels.
    Multi-site entries (e.g. ``T185_Y187``) are expanded so that each
    individual site becomes its own row.  When a site label appears more than
    once (e.g. same gene and site from different peptides), intensities are
    averaged across duplicates.

    The resulting CSV has one row per unique modification site and one column
    per sample, with raw intensity values (no summarisation across samples).
    """
    import pandas as pd

    rows = []
    for col in phospho_df.columns:
        gene = col[0]
        site_str = col[1]
        parsed = parse_phospho_site(site_str)
        if not parsed:
            # Keep raw label when site string is not parseable
            rows.append((f"{gene}_{site_str}", phospho_df[col]))
            continue
        for aa, pos in parsed:
            label = f"{gene}_{aa}{pos}"
            rows.append((label, phospho_df[col]))

    # Build a DataFrame with site labels as columns, then transpose
    labels = [r[0] for r in rows]
    series = [r[1] for r in rows]
    df = pd.DataFrame(dict(enumerate(series)))
    df.columns = labels
    df.index = phospho_df.index

    # Average duplicate site labels (same gene + site from different peptides)
    df = df.T
    df.index.name = "Site"
    df = df.groupby(level=0).mean()

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    df.to_csv(output_path)
    logger.info(
        "Wrote matrix CSV (%d sites x %d samples) to %s",
        len(df),
        len(df.columns),
        output_path,
    )


def write_metadata_csv(clinical_df, cancer_type: str, output_path: str) -> None:
    """Write clinical metadata with an added ``cancer_type`` column."""
    df = clinical_df.copy()
    df["cancer_type"] = cancer_type
    df.index.name = "SampleID"
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    df.to_csv(output_path)
    logger.info("Wrote metadata CSV to %s", output_path)


def melt_proteomics_wide_to_long(proteomics_df):
    """Melt a CPTAC ``get_proteomics()`` wide DataFrame into long format.

    ``get_proteomics()`` returns samples x genes with a 2-level column
    ``MultiIndex`` (``Name``, ``Database_ID``) -- the same column shape
    documented for ``get_phosphoproteomics()`` (see :func:`extract_phospho_sites`),
    minus the ``Site`` level. ``Name`` is the gene symbol; ``Database_ID`` is
    CPTAC's per-column accession, which stays unique even when a protein group
    maps several accessions to the same symbol (isoforms).

    The long format keys each row by ``GeneID`` (= ``Database_ID``, falling back
    to ``Name`` when a column has no accession) rather than by ``Name`` alone,
    so the AnnData builder's ``DataFrame.pivot`` never collides on a duplicate
    ``(SampleID, GeneID)`` pair. ``Name`` is carried through as ``Gene Name``.
    Any residual duplicate ``(SampleID, GeneID)`` rows are collapsed by
    averaging, mirroring :func:`write_matrix_csv`'s duplicate-site handling for
    phospho.

    Returns
    -------
    pd.DataFrame
        Columns ``SampleID``, ``GeneID``, ``Gene Name``, ``Expression`` -- the
        shape ``build_cptac_expression`` reads with its default column names.
    """
    import pandas as pd

    df = proteomics_df
    if isinstance(df.columns, pd.MultiIndex):
        gene_names = list(df.columns.get_level_values(0))
        if "Database_ID" in list(df.columns.names):
            gene_ids = list(df.columns.get_level_values("Database_ID"))
        else:
            gene_ids = list(df.columns.get_level_values(-1))
    else:
        gene_names = list(df.columns)
        gene_ids = list(df.columns)

    # Fall back to the gene symbol when a column has no accession, so a blank
    # Database_ID never collapses two different genes into one GeneID.
    gene_ids = [
        gid if str(gid).strip() and str(gid).lower() != "nan" else name
        for gid, name in zip(gene_ids, gene_names, strict=True)
    ]

    flat = df.copy()
    flat.columns = range(len(flat.columns))
    flat.index.name = "SampleID"
    flat = flat.reset_index()

    long_df = flat.melt(id_vars="SampleID", var_name="_col", value_name="Expression")
    long_df["Gene Name"] = long_df["_col"].map(dict(enumerate(gene_names)))
    long_df["GeneID"] = long_df["_col"].map(dict(enumerate(gene_ids)))
    long_df = long_df.drop(columns="_col")

    long_df = long_df.groupby(["SampleID", "GeneID"], as_index=False).agg(
        {"Gene Name": "first", "Expression": "mean"}
    )
    return long_df[["SampleID", "GeneID", "Gene Name", "Expression"]]


def _tag_tissue_type(clinical_df, sample_ids):
    """Tag each clinical-metadata row with ``tissue_type`` from the sample-id suffix.

    CPTAC encodes normal-tissue samples with a ``.N`` suffix (e.g. ``C3L-00977.N``);
    everything else is tumor. Shared by the phospho and expression downloaders so
    the convention is defined once.
    """
    sample_tissue = {}
    for sid in sample_ids:
        sid_str = str(sid)
        sample_tissue[sid_str] = "normal" if sid_str.endswith(".N") else "tumor"

    clinical_copy = clinical_df.copy()
    clinical_copy["tissue_type"] = clinical_copy.index.map(
        lambda s: sample_tissue.get(str(s), "unknown")
    )
    return clinical_copy


def resolve_raw_dir_inputs(parsed_input, *, cancer_type: Optional[str], prefix: str):
    """Return ``(expression_csv, metadata_csv)`` paths for a per-cancer CPTAC dataset.

    ``parsed_input`` is either a ``{"expression", "metadata"}`` mapping (the direct
    Python API) or a raw download directory -- the ``reprocess`` path for a dataset
    with no ``lifecycle.parse``, which hands the builder the raw dir as-is. In the
    directory case, *cancer_type* selects the ``<prefix>-<cancer_type>-matrix.csv`` /
    ``-metadata.csv`` pair written by that dataset's ``download()``.

    Shared by ``build_cptac_phospho`` (``prefix="cptac-phospho"``) and
    ``build_cptac_expression`` (``prefix="cptac-expression"``).
    """
    if isinstance(parsed_input, Mapping):
        return str(parsed_input["expression"]), str(parsed_input["metadata"])

    raw = str(parsed_input)
    dataset_name = prefix.replace("cptac-", "cptac:", 1)
    if not os.path.isdir(raw):
        raise TypeError(
            f"build_{prefix.replace('-', '_')} expects a "
            "{'expression', 'metadata'} mapping or a raw download directory; "
            f"got {raw!r}"
        )
    if not cancer_type:
        raise ValueError(
            f"{dataset_name} builds one AnnData per cancer type; pass "
            "--plugin-arg cancer_type=<ct> (one of: "
            f"{', '.join(CPTAC_CANCER_TYPES)})"
        )
    expr = os.path.join(raw, f"{prefix}-{cancer_type}-matrix.csv")
    meta = os.path.join(raw, f"{prefix}-{cancer_type}-metadata.csv")
    missing = [p for p in (expr, meta) if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(
            f"Missing {dataset_name} input(s) for cancer_type={cancer_type!r}: "
            f"{missing}. Run `hvantk download {prefix} --cancer-type "
            f"{cancer_type}` (or `hvantk reprocess {dataset_name}`) first."
        )
    return expr, meta


def download_cancer_types(
    dataset_cls,
    raw_dir: str,
    cancer_type: Optional[str] = None,
    overwrite: bool = False,
    *,
    label: str = "CPTAC",
) -> dict:
    """Download one or every CPTAC cancer type via *dataset_cls*, skip-and-continue.

    Shared by the phospho and expression lifecycle ``download_dataset`` entry
    points: run ``dataset_cls(cancer_type=ct).download(raw_dir, overwrite=overwrite)``
    per cancer type, recording any failure instead of raising immediately (``coad``
    and ``ov`` are known upstream failures in ``cptac`` 1.5.14 for phospho; the same
    resilience protects every CPTAC dataset against a bad cancer type). Raises only
    when every cancer type failed.

    Parameters
    ----------
    dataset_cls : type
        ``CPTACPhosphoDataset`` or ``CPTACExpressionDataset``.
    label : str
        Used only in the failure log line and the all-failed error message (e.g.
        ``"CPTAC phospho"``), so each dataset keeps its own wording.
    """
    cancer_types = [cancer_type] if cancer_type else CPTAC_CANCER_TYPES
    results: dict = {}
    failures: dict = {}
    for ct in cancer_types:
        try:
            results[ct] = dataset_cls(cancer_type=ct).download(
                raw_dir, overwrite=overwrite
            )
        except Exception as exc:  # noqa: BLE001 - skip-and-continue per cancer
            logger.warning("%s download failed for %s: %s", label, ct, exc)
            failures[ct] = str(exc)
    if failures:
        results["_failures"] = failures
    succeeded = [k for k in results if k != "_failures"]
    if not succeeded:
        raise RuntimeError(f"All {label} downloads failed: {failures}")
    return results


# ---------------------------------------------------------------------------
# Internal loader
# ---------------------------------------------------------------------------


def _load_cptac_dataset(cancer_type: str):
    """Lazy-import the cptac package and return the dataset object.

    Raises :class:`ImportError` with install instructions when the package
    is not available.
    """
    try:
        import cptac  # noqa: F811
    except ImportError as exc:
        raise ImportError(
            "The 'cptac' package is required for CPTAC data. "
            "Install it with:  pip install cptac"
        ) from exc

    class_name = CPTAC_CANCER_CLASS_MAP[cancer_type]
    dataset_cls = getattr(cptac, class_name)
    return dataset_cls()


# ---------------------------------------------------------------------------
# Dataset class
# ---------------------------------------------------------------------------


@dataclass
class CPTACPhosphoDataset:
    """Dataset class for CPTAC phosphoproteomics data.

    Parameters
    ----------
    cancer_type : str
        One of :data:`~hvantk.skills.cptac.shared.constants.CPTAC_CANCER_TYPES`.
    """

    cancer_type: str

    def __post_init__(self) -> None:
        if self.cancer_type not in CPTAC_CANCER_TYPES:
            raise ValueError(
                f"Unknown cancer type '{self.cancer_type}'. "
                f"Choose from: {', '.join(CPTAC_CANCER_TYPES)}"
            )

    def get_metadata(self) -> dict:
        """Return descriptive metadata about this dataset."""
        return {
            "source": "CPTAC Phosphoproteomics",
            "cancer_type": self.cancer_type,
            "class_name": CPTAC_CANCER_CLASS_MAP[self.cancer_type],
        }

    def download(self, output_dir: str, overwrite: bool = False) -> Dict[str, str]:
        """Download and process CPTAC phospho data into *output_dir*.

        Fetches tumor and normal tissue separately. Produces:
        - ``cptac-phospho-{cancer_type}.tsv``           (combined site table with tissue_type + mean_intensity)
        - ``cptac-phospho-{cancer_type}-tumor.tsv``     (tumor-only sites)
        - ``cptac-phospho-{cancer_type}-normal.tsv``    (normal-only sites, if available)
        - ``cptac-phospho-{cancer_type}-matrix.csv``    (intensity matrix, all samples)
        - ``cptac-phospho-{cancer_type}-metadata.csv``  (clinical metadata with tissue_type)

        Returns
        -------
        dict
            Paths to output files.
        """
        os.makedirs(output_dir, exist_ok=True)
        ct = self.cancer_type

        tsv_path = os.path.join(output_dir, f"cptac-phospho-{ct}.tsv")
        tumor_tsv = os.path.join(output_dir, f"cptac-phospho-{ct}-tumor.tsv")
        normal_tsv = os.path.join(output_dir, f"cptac-phospho-{ct}-normal.tsv")
        matrix_path = os.path.join(output_dir, f"cptac-phospho-{ct}-matrix.csv")
        metadata_path = os.path.join(output_dir, f"cptac-phospho-{ct}-metadata.csv")

        if os.path.exists(tsv_path) and not overwrite:
            logger.info("Output files exist for %s, skipping", ct)
            return {
                "tsv": tsv_path,
                "tumor_tsv": tumor_tsv,
                "normal_tsv": normal_tsv,
                "matrix": matrix_path,
                "metadata": metadata_path,
            }

        logger.info("Loading CPTAC %s dataset...", ct)
        ds = _load_cptac_dataset(ct)

        def _get_phospho(tissue_type="both", source="umich"):
            """Try *source* first; fall back to default if unavailable."""
            try:
                return ds.get_phosphoproteomics(source=source, tissue_type=tissue_type)
            except Exception:
                logger.info(
                    "Source '%s' unavailable for %s, falling back to default",
                    source,
                    ct,
                )
                return ds.get_phosphoproteomics(tissue_type=tissue_type)

        # Fetch tumor phospho data
        logger.info("Fetching tumor phosphoproteomics...")
        tumor_df = _get_phospho(tissue_type="tumor")
        logger.info("Tumor: %d samples x %d columns", *tumor_df.shape)
        tumor_sites = extract_phospho_sites(tumor_df, ct, tissue_type="tumor")
        write_intermediate_tsv(tumor_sites, tumor_tsv)

        # Fetch normal phospho data (may be empty for some cancer types)
        normal_sites = []
        try:
            normal_df = _get_phospho(tissue_type="normal")
            if len(normal_df) > 0:
                logger.info("Normal: %d samples x %d columns", *normal_df.shape)
                normal_sites = extract_phospho_sites(
                    normal_df, ct, tissue_type="normal"
                )
                write_intermediate_tsv(normal_sites, normal_tsv)
            else:
                logger.info("No normal samples available for %s", ct)
        except Exception as e:
            logger.info("No normal tissue data for %s: %s", ct, e)

        # Combined TSV (tumor + normal)
        all_sites = tumor_sites + normal_sites
        write_intermediate_tsv(all_sites, tsv_path)

        # Matrix CSV (all samples — tumor + normal combined)
        both_df = _get_phospho()
        write_matrix_csv(both_df, matrix_path)

        # Clinical metadata with tissue_type per sample
        logger.info("Fetching clinical metadata...")
        try:
            clinical_df = ds.get_clinical(source="mssm")
        except Exception:
            logger.info("Source 'mssm' unavailable for %s clinical, using default", ct)
            clinical_df = ds.get_clinical()

        clinical_copy = _tag_tissue_type(clinical_df, both_df.index)
        write_metadata_csv(clinical_copy, ct, metadata_path)

        n_tumor = len(tumor_sites)
        n_normal = len(normal_sites)
        logger.info(
            "CPTAC %s complete: %d tumor sites, %d normal sites, "
            "%d samples (%d tumor, %d normal)",
            ct,
            n_tumor,
            n_normal,
            len(both_df),
            len(tumor_df),
            len(normal_df) if normal_sites else 0,
        )

        return {
            "tsv": tsv_path,
            "tumor_tsv": tumor_tsv,
            "normal_tsv": normal_tsv,
            "matrix": matrix_path,
            "metadata": metadata_path,
        }


@dataclass
class CPTACExpressionDataset:
    """Dataset class for CPTAC whole-proteome protein-expression data.

    Parameters
    ----------
    cancer_type : str
        One of :data:`~hvantk.skills.cptac.shared.constants.CPTAC_CANCER_TYPES`.
    """

    cancer_type: str

    def __post_init__(self) -> None:
        if self.cancer_type not in CPTAC_CANCER_TYPES:
            raise ValueError(
                f"Unknown cancer type '{self.cancer_type}'. "
                f"Choose from: {', '.join(CPTAC_CANCER_TYPES)}"
            )

    def get_metadata(self) -> dict:
        """Return descriptive metadata about this dataset."""
        return {
            "source": "CPTAC Proteomics",
            "cancer_type": self.cancer_type,
            "class_name": CPTAC_CANCER_CLASS_MAP[self.cancer_type],
        }

    def download(self, output_dir: str, overwrite: bool = False) -> Dict[str, str]:
        """Download and process CPTAC protein-expression data into *output_dir*.

        Fetches whole-proteome abundance (tumor + normal) and clinical metadata.
        Produces:
        - ``cptac-expression-{cancer_type}-matrix.csv``   (long-format table:
          SampleID, GeneID, Gene Name, Expression -- the shape
          ``build_cptac_expression`` reads with its default column names)
        - ``cptac-expression-{cancer_type}-metadata.csv`` (clinical metadata with
          tissue_type)

        Returns
        -------
        dict
            Paths to output files (keys ``matrix``, ``metadata``).
        """
        os.makedirs(output_dir, exist_ok=True)
        ct = self.cancer_type

        matrix_path = os.path.join(output_dir, f"cptac-expression-{ct}-matrix.csv")
        metadata_path = os.path.join(output_dir, f"cptac-expression-{ct}-metadata.csv")

        if (
            os.path.exists(matrix_path)
            and os.path.exists(metadata_path)
            and not overwrite
        ):
            logger.info("Output files exist for %s, skipping", ct)
            return {"matrix": matrix_path, "metadata": metadata_path}

        logger.info("Loading CPTAC %s dataset...", ct)
        ds = _load_cptac_dataset(ct)

        # No fall-back-to-default-source retry here (unlike cptac:phospho's
        # pre-existing _get_phospho/get_clinical fallbacks): cptac 1.5.14's
        # no-source path hits `len(sources_for_data)` on a generator inside
        # cptac/cancers/cancer.py get_dataframe() (~line 682) and raises an
        # opaque `TypeError: object of type 'generator' has no len()` instead
        # of the real "no source available" condition. Always pass an
        # explicit source and, on failure, raise a clear error naming the
        # cancer type, data type, and source instead of retrying into that bug.
        proteomics_source = "umich"
        logger.info("Fetching proteomics (source=%s)...", proteomics_source)
        try:
            proteomics_df = ds.get_proteomics(
                source=proteomics_source, tissue_type="both"
            )
        except Exception as exc:
            raise RuntimeError(
                f"CPTAC proteomics fetch failed for cancer_type={ct!r}, "
                f"data_type='proteomics', source={proteomics_source!r}: {exc}"
            ) from exc
        logger.info("Proteomics: %d samples x %d columns", *proteomics_df.shape)

        clinical_source = "mssm"
        logger.info("Fetching clinical metadata (source=%s)...", clinical_source)
        try:
            clinical_df = ds.get_clinical(source=clinical_source)
        except Exception as exc:
            raise RuntimeError(
                f"CPTAC clinical-metadata fetch failed for cancer_type={ct!r}, "
                f"data_type='clinical', source={clinical_source!r}: {exc}"
            ) from exc

        # Both fetches succeeded; compute both outputs before writing anything.
        long_df = melt_proteomics_wide_to_long(proteomics_df)
        clinical_copy = _tag_tissue_type(clinical_df, proteomics_df.index)

        # Write to temp files first and publish (os.replace, atomic within one
        # filesystem) only once BOTH writes have succeeded. Publishing the matrix
        # as soon as it was ready -- the previous behavior -- could leave a fresh
        # matrix beside a stale (or missing) metadata file if the clinical fetch or
        # write then failed, and a bare `df.to_csv(final_path)` is not atomic either
        # (a crash or error mid-write leaves a truncated file under the final name).
        matrix_fd, matrix_tmp = tempfile.mkstemp(
            dir=output_dir, prefix=f".{os.path.basename(matrix_path)}."
        )
        os.close(matrix_fd)
        metadata_fd, metadata_tmp = tempfile.mkstemp(
            dir=output_dir, prefix=f".{os.path.basename(metadata_path)}."
        )
        os.close(metadata_fd)
        try:
            long_df.to_csv(matrix_tmp, index=False)
            write_metadata_csv(clinical_copy, ct, metadata_tmp)
            # mkstemp always creates mode 0600; match what a normal file
            # create (or the existing file, on overwrite) would get, so a
            # shared raw dir stays group-readable. Mirrors
            # drift_runner.write_fingerprint's atomic-publish idiom.
            for tmp, final in (
                (matrix_tmp, matrix_path),
                (metadata_tmp, metadata_path),
            ):
                if os.path.exists(final):
                    shutil.copymode(final, tmp)
                else:
                    os.chmod(tmp, 0o666 & ~_current_umask())
            os.replace(matrix_tmp, matrix_path)
            os.replace(metadata_tmp, metadata_path)
        finally:
            for tmp in (matrix_tmp, metadata_tmp):
                try:
                    os.unlink(tmp)
                except FileNotFoundError:
                    pass  # already published (renamed) or never written

        logger.info(
            "Wrote expression matrix (%d rows) to %s", len(long_df), matrix_path
        )
        logger.info(
            "CPTAC %s expression complete: %d samples, %d (sample, gene) rows",
            ct,
            len(proteomics_df),
            len(long_df),
        )

        return {"matrix": matrix_path, "metadata": metadata_path}


def parse_raw_dir(
    raw_dir: str, output_path: str, cancer_type: Optional[str] = None, **kwargs
) -> str:
    """Lifecycle ``parse`` entry point for the plugin loader.

    Per the ``DatasetSpec`` contract in :mod:`hvantk.core.plugin_api`, a
    lifecycle ``parse_fn`` reads files written by ``download_fn`` from
    ``raw_dir`` and writes an intermediate representation to ``output_path``.

    CPTAC data lives in pre-parsed per-cancer-type TSVs already written by
    ``CPTACPhosphoDataset.download``. This function consolidates any
    ``cptac-phospho-*.tsv`` files found under ``raw_dir`` into a single
    pan-cancer TSV at ``output_path``.

    Parameters
    ----------
    raw_dir : str
        Directory containing one or more ``cptac-phospho-<cancer_type>.tsv``
        files produced by :class:`CPTACPhosphoDataset`.
    output_path : str
        Destination consolidated TSV path.
    cancer_type : str, optional
        If provided, only include the matching cancer type's TSV (this is the
        single-source case where ``parse_raw_dir`` simply copies the existing
        per-cancer TSV to ``output_path``).
    **kwargs
        Reserved for future lifecycle keyword arguments.

    Returns
    -------
    str
        The output TSV path.
    """
    if cancer_type is not None:
        pattern_names = [f"cptac-phospho-{cancer_type}.tsv"]
    else:
        pattern_names = [
            name
            for name in os.listdir(raw_dir)
            if name.startswith("cptac-phospho-")
            and name.endswith(".tsv")
            and "-tumor" not in name
            and "-normal" not in name
            and not name.endswith("-pancancer.tsv")
        ]

    if not pattern_names:
        raise FileNotFoundError(
            f"No cptac-phospho-*.tsv files found in {raw_dir!r}. "
            "Run `hvantk download cptac-phospho` first."
        )

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", newline="") as fout:
        writer = csv.DictWriter(
            fout,
            fieldnames=_TSV_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for name in pattern_names:
            tsv_path = os.path.join(raw_dir, name)
            with open(tsv_path) as fin:
                reader = csv.DictReader(fin, delimiter="\t")
                for row in reader:
                    writer.writerow(row)

    logger.info(
        "Consolidated %d cancer-type TSV(s) into %s", len(pattern_names), output_path
    )
    return output_path
