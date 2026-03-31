"""
CPTAC Phospho dataset handling.

This module provides functions and a dataset class for downloading and parsing
phosphoproteomics data from CPTAC (Clinical Proteomic Tumor Analysis
Consortium) via the ``cptac`` Python package.

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
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from hvantk.ptm.constants import CPTAC_CANCER_TYPES, CPTAC_CANCER_CLASS_MAP

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

    The multi-index columns are flattened to ``Gene_Site`` labels and the
    DataFrame is transposed so that sites become rows and samples become
    columns.
    """
    df = phospho_df.copy()
    df.columns = [f"{c[0]}_{c[1]}" for c in df.columns]
    df = df.T
    df.index.name = "Site"
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    df.to_csv(output_path)
    logger.info("Wrote matrix CSV to %s", output_path)


def write_metadata_csv(clinical_df, cancer_type: str, output_path: str) -> None:
    """Write clinical metadata with an added ``cancer_type`` column."""
    df = clinical_df.copy()
    df["cancer_type"] = cancer_type
    df.index.name = "SampleID"
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    df.to_csv(output_path)
    logger.info("Wrote metadata CSV to %s", output_path)


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
        One of :data:`~hvantk.ptm.constants.CPTAC_CANCER_TYPES`.
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
                "tsv": tsv_path, "tumor_tsv": tumor_tsv,
                "normal_tsv": normal_tsv, "matrix": matrix_path,
                "metadata": metadata_path,
            }

        logger.info("Loading CPTAC %s dataset...", ct)
        ds = _load_cptac_dataset(ct)

        # Fetch tumor phospho data
        logger.info("Fetching tumor phosphoproteomics...")
        tumor_df = ds.get_phosphoproteomics(source="umich", tissue_type="tumor")
        logger.info("Tumor: %d samples x %d columns", *tumor_df.shape)
        tumor_sites = extract_phospho_sites(tumor_df, ct, tissue_type="tumor")
        write_intermediate_tsv(tumor_sites, tumor_tsv)

        # Fetch normal phospho data (may be empty for some cancer types)
        normal_sites = []
        try:
            normal_df = ds.get_phosphoproteomics(source="umich", tissue_type="normal")
            if len(normal_df) > 0:
                logger.info("Normal: %d samples x %d columns", *normal_df.shape)
                normal_sites = extract_phospho_sites(normal_df, ct, tissue_type="normal")
                write_intermediate_tsv(normal_sites, normal_tsv)
            else:
                logger.info("No normal samples available for %s", ct)
        except Exception as e:
            logger.info("No normal tissue data for %s: %s", ct, e)

        # Combined TSV (tumor + normal)
        all_sites = tumor_sites + normal_sites
        write_intermediate_tsv(all_sites, tsv_path)

        # Matrix CSV (all samples — tumor + normal combined)
        both_df = ds.get_phosphoproteomics(source="umich")
        write_matrix_csv(both_df, matrix_path)

        # Clinical metadata with tissue_type per sample
        logger.info("Fetching clinical metadata...")
        clinical_df = ds.get_clinical(source="mssm")

        # Tag samples with tissue_type based on sample ID suffix
        sample_tissue = {}
        for sid in both_df.index:
            sid_str = str(sid)
            if sid_str.endswith(".N"):
                sample_tissue[sid_str] = "normal"
            else:
                sample_tissue[sid_str] = "tumor"

        clinical_copy = clinical_df.copy()
        clinical_copy["tissue_type"] = clinical_copy.index.map(
            lambda s: sample_tissue.get(str(s), "unknown")
        )
        write_metadata_csv(clinical_copy, ct, metadata_path)

        n_tumor = len(tumor_sites)
        n_normal = len(normal_sites)
        logger.info(
            "CPTAC %s complete: %d tumor sites, %d normal sites, "
            "%d samples (%d tumor, %d normal)",
            ct, n_tumor, n_normal, len(both_df),
            len(tumor_df), len(normal_df) if normal_sites else 0,
        )

        return {
            "tsv": tsv_path, "tumor_tsv": tumor_tsv,
            "normal_tsv": normal_tsv, "matrix": matrix_path,
            "metadata": metadata_path,
        }
