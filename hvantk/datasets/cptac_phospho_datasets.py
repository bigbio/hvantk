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
) -> List[dict]:
    """Extract phospho sites from a CPTAC multi-index DataFrame.

    Parameters
    ----------
    phospho_df : pd.DataFrame
        DataFrame with a :class:`~pandas.MultiIndex` on columns whose first
        two levels are *Gene* and *Site*.
    cancer_type : str
        CPTAC cancer type identifier (e.g. ``"brca"``).

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
                # Merge observations from multiple columns mapping to the
                # same site (unlikely but defensive).
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
                    "mean_intensity": mean_val,
                    "source_db": "CPTAC",
                    "evidence_type": "mass_spectrometry",
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

        Produces three files:
        - ``cptac-phospho-{cancer_type}.tsv``  (intermediate site table)
        - ``cptac-phospho-{cancer_type}-matrix.csv``   (intensity matrix)
        - ``cptac-phospho-{cancer_type}-metadata.csv`` (clinical metadata)

        Returns
        -------
        dict
            Paths to output files: {"tsv": ..., "matrix": ..., "metadata": ...}
        """
        os.makedirs(output_dir, exist_ok=True)

        tsv_path = os.path.join(output_dir, f"cptac-phospho-{self.cancer_type}.tsv")
        matrix_path = os.path.join(output_dir, f"cptac-phospho-{self.cancer_type}-matrix.csv")
        metadata_path = os.path.join(output_dir, f"cptac-phospho-{self.cancer_type}-metadata.csv")

        if all(os.path.exists(p) for p in [tsv_path, matrix_path, metadata_path]) and not overwrite:
            logger.info("All output files exist for %s, skipping", self.cancer_type)
            return {"tsv": tsv_path, "matrix": matrix_path, "metadata": metadata_path}

        logger.info("Loading CPTAC %s dataset...", self.cancer_type)
        ds = _load_cptac_dataset(self.cancer_type)

        logger.info("Fetching phosphoproteomics data...")
        # Specify source='umich' to avoid a bug in cptac when multiple
        # sources exist and no source is specified (generator has no len()).
        phospho_df = ds.get_phosphoproteomics(source="umich")
        logger.info("Phospho DataFrame: %d samples x %d columns", *phospho_df.shape)

        logger.info("Fetching clinical metadata...")
        # Specify source='mssm' to avoid the same multi-source bug.
        clinical_df = ds.get_clinical(source="mssm")

        sites = extract_phospho_sites(phospho_df, self.cancer_type)
        write_intermediate_tsv(sites, tsv_path)
        write_matrix_csv(phospho_df, matrix_path)
        write_metadata_csv(clinical_df, self.cancer_type, metadata_path)

        logger.info(
            "CPTAC %s download complete: %d sites from %d samples",
            self.cancer_type,
            len(sites),
            len(phospho_df),
        )
        return {"tsv": tsv_path, "matrix": matrix_path, "metadata": metadata_path}
