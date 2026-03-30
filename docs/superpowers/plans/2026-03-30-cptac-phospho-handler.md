# CPTAC Phospho Handler Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add CPTAC phosphoproteomics as a data source in the PTM pipeline, using the `cptac` Python package to extract phospho sites with per-sample intensities from 7+ cancer types, producing both a site-level summary TSV for genomic coordinate mapping and a phospho MatrixTable (sites × samples).

**Architecture:** New dataset class `CPTACPhosphoDataset` wraps the `cptac` Python package, parses phospho site names from multi-index DataFrames, aggregates per-site stats, and writes both an intermediate TSV (for the PTM pipeline) and a matrix CSV + metadata CSV (for the MatrixTable builder). Reuses the existing GTF mapper via `gene_mane` fallback and follows the PeptideAtlas concatenation pattern for pipeline integration.

**Tech Stack:** Python 3.10+, `cptac` package (optional dep, lazy-imported), Click (CLI), pandas, csv, pytest, Hail (MatrixTable builder)

**Spec:** `docs/superpowers/specs/2026-03-30-cptac-phospho-handler-design.md`

---

## File Structure

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `hvantk/datasets/cptac_phospho_datasets.py` | Dataset class: wrap cptac package, parse sites, aggregate, write TSV + matrix + metadata |
| Create | `hvantk/commands/cptac_phospho_downloader.py` | CLI `hvantk download cptac-phospho` |
| Create | `hvantk/tests/test_cptac_phospho.py` | Unit tests with mock DataFrames (no real cptac dependency) |
| Modify | `hvantk/ptm/constants.py` | Add CPTAC constants (cancer types, phospho AA map) |
| Modify | `hvantk/ptm/pipeline.py` | Add `cptac_tsv` to PTMBuildConfig, concatenation in ptm_build_pipeline() |
| Modify | `hvantk/commands/ptm_cli.py` | Add `--cptac-tsv` option to ptm build |
| Modify | `hvantk/commands/download_cli.py` | Register cptac-phospho downloader |
| Modify | `hvantk/tables/matrix_builders.py` | Add `build_cptac_phospho_mt()` |
| Modify | `hvantk/tables/registry.py` | Register "cptac-phospho" in MATRIX_BUILDERS |
| Modify | `hvantk/commands/make_matrix_cli.py` | Add `hvantk mkmatrix cptac-phospho` command |

---

### Task 1: Add CPTAC Constants

**Files:**
- Modify: `hvantk/ptm/constants.py`

- [ ] **Step 1: Add CPTAC constants after the PeptideAtlas block**

In `hvantk/ptm/constants.py`, after the PeptideAtlas constants (around line 20), add:

```python
# CPTAC Phospho (via cptac Python package)
CPTAC_CANCER_TYPES = [
    "brca", "ccrcc", "colon", "endometrial", "gbm",
    "hnscc", "lscc", "luad", "ov", "pdac", "ucec",
]
CPTAC_CANCER_CLASS_MAP = {
    "brca": "Brca",
    "ccrcc": "Ccrcc",
    "colon": "Colon",
    "endometrial": "Endometrial",
    "gbm": "Gbm",
    "hnscc": "Hnscc",
    "lscc": "Lscc",
    "luad": "Luad",
    "ov": "Ov",
    "pdac": "Pdac",
    "ucec": "Ucec",
}
```

- [ ] **Step 2: Commit**

```bash
git add hvantk/ptm/constants.py
git commit -m "feat: add CPTAC cancer type constants"
```

---

### Task 2: Create CPTAC Phospho Dataset Class — Tests First

**Files:**
- Create: `hvantk/tests/test_cptac_phospho.py`
- Create: `hvantk/datasets/cptac_phospho_datasets.py`

- [ ] **Step 1: Write failing tests**

Create `hvantk/tests/test_cptac_phospho.py`:

```python
"""Tests for hvantk.datasets.cptac_phospho_datasets module.

Covers:
1. Site name parsing (single site, multi-site, edge cases)
2. Phospho site extraction from mock DataFrame
3. Intermediate TSV output format
4. Matrix CSV and metadata CSV output
5. Dataset class construction
"""

import csv
import os

import numpy as np
import pandas as pd
import pytest


# ---------- Test 1: Site name parsing ----------


@pytest.mark.parametrize(
    "site_str, expected",
    [
        ("S65", [("S", 65)]),
        ("T185", [("T", 185)]),
        ("Y243", [("Y", 243)]),
        ("T185_Y187", [("T", 185), ("Y", 187)]),
        ("S15_S20_T25", [("S", 15), ("S", 20), ("T", 25)]),
        ("", []),
        ("X100", []),  # Non-S/T/Y
    ],
)
def test_parse_phospho_site(site_str, expected):
    """Site name parser handles single, multi, and edge cases."""
    from hvantk.datasets.cptac_phospho_datasets import parse_phospho_site

    result = parse_phospho_site(site_str)
    assert result == expected


# ---------- Fixtures ----------


@pytest.fixture
def mock_phospho_df():
    """Mock phosphoproteomics DataFrame mimicking cptac package output.

    Multi-index columns: (Gene, Site, Peptide, Database_ID)
    Rows: samples
    Values: log2 intensities (some NaN)
    """
    columns = pd.MultiIndex.from_tuples(
        [
            ("TP53", "S315", "SPQPKKKPLDGEpS", "NP_000537.3"),
            ("TP53", "S6", "MEEPQpSDPSVEPPL", "NP_000537.3"),
            ("MAPK1", "T185_Y187", "VADPDHDHTGFLpTEpYVATR", "NP_002736.3"),
            ("EIF4EBP1", "S65", "RVpSGGEELGS", "NP_004086.1"),
        ],
        names=["Gene", "Site", "Peptide", "Database_ID"],
    )
    data = np.array(
        [
            [1.5, 2.0, 0.5, np.nan],  # Sample_01
            [np.nan, 1.8, -0.3, 3.2],  # Sample_02
            [0.7, np.nan, 1.1, 2.5],  # Sample_03
        ]
    )
    return pd.DataFrame(data, index=["Sample_01", "Sample_02", "Sample_03"], columns=columns)


@pytest.fixture
def mock_clinical_df():
    """Mock clinical DataFrame."""
    return pd.DataFrame(
        {
            "Sample_type": ["Tumor", "Tumor", "Normal"],
            "Age": [55, 63, 48],
            "Gender": ["Female", "Male", "Female"],
        },
        index=["Sample_01", "Sample_02", "Sample_03"],
    )


# ---------- Test 2: Phospho site extraction ----------


def test_extract_phospho_sites(mock_phospho_df):
    """Phospho sites extracted with correct aggregation."""
    from hvantk.datasets.cptac_phospho_datasets import extract_phospho_sites

    sites = extract_phospho_sites(mock_phospho_df, cancer_type="brca")

    # TP53 S315: 2 non-NaN values (Sample_01=1.5, Sample_03=0.7)
    s315 = next(s for s in sites if s["gene_symbol"] == "TP53" and s["position"] == 315)
    assert s315["amino_acid"] == "S"
    assert s315["n_observations"] == 2
    assert abs(s315["mean_intensity"] - (1.5 + 0.7) / 2) < 0.01
    assert s315["cancer_type"] == "brca"

    # MAPK1 T185_Y187 → two separate sites (T185 and Y187)
    mapk_sites = [s for s in sites if s["gene_symbol"] == "MAPK1"]
    positions = {s["position"] for s in mapk_sites}
    assert positions == {185, 187}

    # All sites should have source_db = CPTAC
    assert all(s["source_db"] == "CPTAC" for s in sites)


# ---------- Test 3: Intermediate TSV output ----------


def test_write_intermediate_tsv(mock_phospho_df, tmp_path):
    """Intermediate TSV has correct columns and values."""
    from hvantk.datasets.cptac_phospho_datasets import (
        extract_phospho_sites,
        write_intermediate_tsv,
    )

    sites = extract_phospho_sites(mock_phospho_df, cancer_type="brca")
    output_path = str(tmp_path / "cptac_phospho.tsv")
    write_intermediate_tsv(sites, output_path)

    with open(output_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    # 5 sites: TP53 S315, TP53 S6, MAPK1 T185, MAPK1 Y187, EIF4EBP1 S65
    assert len(rows) == 5

    expected_cols = {
        "accession", "gene_symbol", "position", "description",
        "amino_acid", "ensembl_xrefs", "sequence_length",
        "n_observations", "source_db", "evidence_type",
    }
    assert set(reader.fieldnames) == expected_cols

    row_s315 = next(r for r in rows if r["gene_symbol"] == "TP53" and r["position"] == "315")
    assert row_s315["source_db"] == "CPTAC"
    assert row_s315["evidence_type"] == "mass_spectrometry"
    assert row_s315["description"] == "Phosphoserine"


# ---------- Test 4: Matrix CSV output ----------


def test_write_matrix_csv(mock_phospho_df, tmp_path):
    """Matrix CSV has sites as rows and samples as columns."""
    from hvantk.datasets.cptac_phospho_datasets import write_matrix_csv

    output_path = str(tmp_path / "matrix.csv")
    write_matrix_csv(mock_phospho_df, output_path)

    df = pd.read_csv(output_path, index_col=0)
    # Columns should be sample IDs
    assert set(df.columns) == {"Sample_01", "Sample_02", "Sample_03"}
    # Rows should be site identifiers (Gene_Site)
    assert "TP53_S315" in df.index.tolist()


# ---------- Test 5: Metadata CSV output ----------


def test_write_metadata_csv(mock_clinical_df, tmp_path):
    """Metadata CSV preserves clinical columns."""
    from hvantk.datasets.cptac_phospho_datasets import write_metadata_csv

    output_path = str(tmp_path / "metadata.csv")
    write_metadata_csv(mock_clinical_df, "brca", output_path)

    df = pd.read_csv(output_path, index_col=0)
    assert "cancer_type" in df.columns
    assert "Sample_type" in df.columns
    assert df.loc["Sample_01", "cancer_type"] == "brca"


# ---------- Test 6: Dataset class ----------


def test_dataset_class():
    """CPTACPhosphoDataset construction and metadata."""
    from hvantk.datasets.cptac_phospho_datasets import CPTACPhosphoDataset

    dataset = CPTACPhosphoDataset(cancer_type="brca")
    assert dataset.cancer_type == "brca"
    meta = dataset.get_metadata()
    assert meta["source"] == "CPTAC Phosphoproteomics"
    assert meta["cancer_type"] == "brca"


def test_dataset_invalid_cancer():
    """Invalid cancer type raises ValueError."""
    from hvantk.datasets.cptac_phospho_datasets import CPTACPhosphoDataset

    with pytest.raises(ValueError, match="Unknown cancer type"):
        CPTACPhosphoDataset(cancer_type="invalid_cancer")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest hvantk/tests/test_cptac_phospho.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement the dataset class**

Create `hvantk/datasets/cptac_phospho_datasets.py`:

```python
"""
CPTAC Phosphoproteomics dataset handling.

Wraps the ``cptac`` Python package to extract phosphorylation site
positions with per-sample log2 intensities from CPTAC cancer cohorts.
Produces an intermediate TSV for PTM pipeline integration and a
matrix CSV + metadata CSV for MatrixTable construction.

Example usage:
    dataset = CPTACPhosphoDataset(cancer_type="brca")
    dataset.download("/data/cptac/")

    # Or work with a pre-loaded DataFrame:
    sites = extract_phospho_sites(phospho_df, cancer_type="brca")
    write_intermediate_tsv(sites, "cptac_phospho.tsv")
"""

import csv
import logging
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from hvantk.ptm.constants import CPTAC_CANCER_TYPES, CPTAC_CANCER_CLASS_MAP

logger = logging.getLogger(__name__)

# Phospho amino acid descriptions (shared with PeptideAtlas handler)
_PHOSPHO_AA_DESC = {
    "S": "Phosphoserine",
    "T": "Phosphothreonine",
    "Y": "Phosphotyrosine",
}

# Intermediate TSV columns
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

# Pattern to parse site names like "S65", "T185", "Y243"
_SITE_RE = re.compile(r"([STY])(\d+)")


def parse_phospho_site(site_str: str) -> List[Tuple[str, int]]:
    """Parse a CPTAC phospho site string into (amino_acid, position) tuples.

    Single site: ``"S65"`` → ``[("S", 65)]``
    Multi-site: ``"T185_Y187"`` → ``[("T", 185), ("Y", 187)]``

    Only S, T, Y residues are accepted.
    """
    if not site_str:
        return []
    return [(m.group(1), int(m.group(2))) for m in _SITE_RE.finditer(site_str)]


def extract_phospho_sites(
    phospho_df: pd.DataFrame,
    cancer_type: str,
) -> List[dict]:
    """Extract phospho sites from a cptac phosphoproteomics DataFrame.

    Parameters
    ----------
    phospho_df : pd.DataFrame
        DataFrame from ``dataset.get_phosphoproteomics()``.
        Multi-index columns: (Gene, Site, Peptide, Database_ID).
        Rows: samples. Values: log2 intensities.
    cancer_type : str
        Cancer type label (e.g., "brca").

    Returns
    -------
    list of dict
        Each dict has keys: gene_symbol, position, amino_acid, description,
        n_observations, mean_intensity, cancer_type, source_db, evidence_type,
        accession, ensembl_xrefs, sequence_length.
    """
    sites = []
    seen = set()

    for col_idx in range(len(phospho_df.columns)):
        col = phospho_df.columns[col_idx]
        # Multi-index: (Gene, Site, Peptide, Database_ID)
        if isinstance(col, tuple) and len(col) >= 2:
            gene = col[0]
            site_str = col[1]
        else:
            continue

        values = phospho_df.iloc[:, col_idx]
        non_nan = values.dropna()
        n_obs = len(non_nan)
        mean_int = float(non_nan.mean()) if n_obs > 0 else 0.0

        parsed = parse_phospho_site(site_str)
        for aa, pos in parsed:
            key = (gene, aa, pos)
            if key in seen:
                continue
            seen.add(key)

            sites.append(
                {
                    "gene_symbol": gene,
                    "position": pos,
                    "amino_acid": aa,
                    "description": _PHOSPHO_AA_DESC.get(aa, "Phosphorylation"),
                    "n_observations": n_obs,
                    "mean_intensity": mean_int,
                    "cancer_type": cancer_type,
                    "source_db": "CPTAC",
                    "evidence_type": "mass_spectrometry",
                    "accession": "",
                    "ensembl_xrefs": "",
                    "sequence_length": "",
                }
            )

    logger.info("Extracted %d phospho sites from %s", len(sites), cancer_type)
    return sites


def write_intermediate_tsv(sites: List[dict], output_path: str) -> str:
    """Write phospho sites to an intermediate TSV for the PTM pipeline.

    Parameters
    ----------
    sites : list of dict
        Output from :func:`extract_phospho_sites`.
    output_path : str
        Destination file path.

    Returns
    -------
    str
        The output file path.
    """
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=_TSV_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for site in sites:
            writer.writerow(site)

    logger.info("Wrote %d phospho sites to %s", len(sites), output_path)
    return output_path


def write_matrix_csv(
    phospho_df: pd.DataFrame,
    output_path: str,
) -> str:
    """Write phospho intensities as a sites-by-samples matrix CSV.

    Rows are site identifiers (Gene_Site), columns are sample IDs,
    values are log2 intensities.

    Parameters
    ----------
    phospho_df : pd.DataFrame
        DataFrame from ``dataset.get_phosphoproteomics()``.
    output_path : str
        Destination file path.

    Returns
    -------
    str
        The output file path.
    """
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    # Transpose: samples become columns, sites become rows
    df_t = phospho_df.T.copy()

    # Flatten multi-index to "Gene_Site" row labels
    if isinstance(df_t.index, pd.MultiIndex):
        df_t.index = [f"{gene}_{site}" for gene, site, *_ in df_t.index]
    df_t.index.name = "SiteID"

    df_t.to_csv(output_path)
    logger.info("Wrote matrix CSV (%d sites x %d samples) to %s",
                len(df_t), len(df_t.columns), output_path)
    return output_path


def write_metadata_csv(
    clinical_df: pd.DataFrame,
    cancer_type: str,
    output_path: str,
) -> str:
    """Write sample clinical metadata as CSV.

    Adds a ``cancer_type`` column and ensures the index is named ``SampleID``.

    Parameters
    ----------
    clinical_df : pd.DataFrame
        DataFrame from ``dataset.get_clinical()``.
    cancer_type : str
        Cancer type label.
    output_path : str
        Destination file path.

    Returns
    -------
    str
        The output file path.
    """
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    df = clinical_df.copy()
    df["cancer_type"] = cancer_type
    df.index.name = "SampleID"
    df.to_csv(output_path)

    logger.info("Wrote metadata for %d samples to %s", len(df), output_path)
    return output_path


def _load_cptac_dataset(cancer_type: str):
    """Lazy-import cptac and load a cancer dataset.

    Raises ImportError with install instructions if cptac is not installed.
    """
    try:
        import cptac
    except ImportError:
        raise ImportError(
            "The 'cptac' package is required for CPTAC data access. "
            "Install it with: pip install cptac"
        ) from None

    class_name = CPTAC_CANCER_CLASS_MAP.get(cancer_type)
    if class_name is None:
        raise ValueError(
            f"Unknown cancer type: {cancer_type}. "
            f"Available: {', '.join(CPTAC_CANCER_TYPES)}"
        )

    dataset_cls = getattr(cptac, class_name)
    return dataset_cls()


@dataclass
class CPTACPhosphoDataset:
    """Represents a CPTAC phosphoproteomics dataset for one cancer type.

    Attributes:
        cancer_type: Cancer type identifier (e.g., "brca").
    """

    cancer_type: str

    def __post_init__(self):
        if self.cancer_type not in CPTAC_CANCER_TYPES:
            raise ValueError(
                f"Unknown cancer type: {self.cancer_type}. "
                f"Available: {', '.join(CPTAC_CANCER_TYPES)}"
            )

    def download(self, output_dir: str, overwrite: bool = False) -> Dict[str, str]:
        """Download and process CPTAC phospho data for this cancer type.

        Produces three output files:
        - Intermediate TSV for PTM pipeline
        - Matrix CSV (sites × samples)
        - Metadata CSV (sample clinical info)

        Parameters
        ----------
        output_dir : str
            Directory to save output files.
        overwrite : bool
            If True, overwrite existing files.

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
        dataset = _load_cptac_dataset(self.cancer_type)

        logger.info("Fetching phosphoproteomics data...")
        phospho_df = dataset.get_phosphoproteomics()
        logger.info("Phospho DataFrame: %d samples x %d columns", *phospho_df.shape)

        logger.info("Fetching clinical metadata...")
        clinical_df = dataset.get_clinical()

        # Write outputs
        sites = extract_phospho_sites(phospho_df, self.cancer_type)
        write_intermediate_tsv(sites, tsv_path)
        write_matrix_csv(phospho_df, matrix_path)
        write_metadata_csv(clinical_df, self.cancer_type, metadata_path)

        return {"tsv": tsv_path, "matrix": matrix_path, "metadata": metadata_path}

    def get_metadata(self) -> dict:
        """Return metadata about this dataset."""
        return {
            "source": "CPTAC Phosphoproteomics",
            "cancer_type": self.cancer_type,
            "description": (
                f"Phosphoproteomics data from the CPTAC {self.cancer_type.upper()} cohort. "
                "Sites are log2 intensity ratios from TMT/iTRAQ labeling."
            ),
        }

    def __str__(self) -> str:
        return f"CPTACPhosphoDataset(cancer_type={self.cancer_type})"
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest hvantk/tests/test_cptac_phospho.py -v`
Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add hvantk/datasets/cptac_phospho_datasets.py hvantk/tests/test_cptac_phospho.py
git commit -m "feat: add CPTAC phospho dataset class with site parsing and tests"
```

---

### Task 3: Create CLI Downloader and Register It

**Files:**
- Create: `hvantk/commands/cptac_phospho_downloader.py`
- Modify: `hvantk/commands/download_cli.py`

- [ ] **Step 1: Create the downloader CLI command**

Create `hvantk/commands/cptac_phospho_downloader.py`:

```python
"""
CLI command to download CPTAC phosphoproteomics data.

Example:
    hvantk download cptac-phospho --cancer-type brca -o data/cptac/
    hvantk download cptac-phospho --all -o data/cptac/
    hvantk download cptac-phospho --list-cancers
"""

import logging

import click

from hvantk.ptm.constants import CPTAC_CANCER_TYPES

logger = logging.getLogger(__name__)


@click.command(name="cptac-phospho")
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    default=None,
    help="Directory to save CPTAC phospho output files",
)
@click.option(
    "--cancer-type",
    type=click.Choice(CPTAC_CANCER_TYPES, case_sensitive=False),
    default=None,
    help="Cancer type to download",
)
@click.option(
    "--all",
    "download_all",
    is_flag=True,
    help="Download all cancer types and produce pan-cancer output",
)
@click.option(
    "--list-cancers",
    is_flag=True,
    help="List available cancer types and exit",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing files",
)
@click.pass_context
def cptac_phospho_downloader(ctx, output_dir, cancer_type, download_all, list_cancers, overwrite):
    """Download CPTAC phosphoproteomics data.

    Downloads phospho site intensities from the CPTAC Python package,
    extracts site positions, and writes files for PTM pipeline and
    MatrixTable construction.

    \b
    Requires: pip install cptac
    """
    if list_cancers:
        click.echo("Available CPTAC cancer types:")
        for ct in CPTAC_CANCER_TYPES:
            click.echo(f"  {ct}")
        return

    if not output_dir:
        click.echo("Error: --output-dir is required (unless using --list-cancers)", err=True)
        ctx.exit(1)

    if not cancer_type and not download_all:
        click.echo(
            "Error: specify --cancer-type or --all. "
            "Use --list-cancers to see available types.",
            err=True,
        )
        ctx.exit(1)

    try:
        from hvantk.datasets.cptac_phospho_datasets import CPTACPhosphoDataset

        cancer_types = CPTAC_CANCER_TYPES if download_all else [cancer_type]

        all_tsv_paths = []
        for ct in cancer_types:
            click.echo(f"Processing {ct}...")
            dataset = CPTACPhosphoDataset(cancer_type=ct)
            paths = dataset.download(output_dir, overwrite=overwrite)
            click.echo(f"  TSV: {paths['tsv']}")
            click.echo(f"  Matrix: {paths['matrix']}")
            click.echo(f"  Metadata: {paths['metadata']}")
            all_tsv_paths.append(paths["tsv"])

        # Pan-cancer merge if --all
        if download_all and len(all_tsv_paths) > 1:
            import csv

            pancancer_path = os.path.join(output_dir, "cptac-phospho-pancancer.tsv")
            click.echo(f"Merging {len(all_tsv_paths)} cancer types into {pancancer_path}...")

            from hvantk.datasets.cptac_phospho_datasets import _TSV_COLUMNS

            with open(pancancer_path, "w", newline="") as fout:
                writer = csv.DictWriter(
                    fout, fieldnames=_TSV_COLUMNS, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                for tsv_path in all_tsv_paths:
                    with open(tsv_path) as fin:
                        reader = csv.DictReader(fin, delimiter="\t")
                        for row in reader:
                            writer.writerow(row)

            click.echo(f"Pan-cancer TSV: {pancancer_path}")

    except ImportError as e:
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
    except Exception as e:
        logger.exception(f"Download failed: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
```

- [ ] **Step 2: Register the downloader in download_cli.py**

In `hvantk/commands/download_cli.py`, add after the peptideatlas import (line 11):

```python
from hvantk.commands.cptac_phospho_downloader import cptac_phospho_downloader
```

Add after the peptideatlas registration (line 26):

```python
download_group.add_command(cptac_phospho_downloader, "cptac-phospho")
```

- [ ] **Step 3: Add missing `import os` to the downloader**

The `os.path.join` call in the pan-cancer merge needs `import os`. Add at the top of `cptac_phospho_downloader.py`:

```python
import os
```

- [ ] **Step 4: Verify CLI registration**

Run: `python -c "from hvantk.commands.download_cli import download_group; from click.testing import CliRunner; r = CliRunner().invoke(download_group, ['cptac-phospho', '--list-cancers']); print(r.output); assert 'brca' in r.output"`

- [ ] **Step 5: Commit**

```bash
git add hvantk/commands/cptac_phospho_downloader.py hvantk/commands/download_cli.py
git commit -m "feat: add cptac-phospho download CLI command"
```

---

### Task 4: Integrate CPTAC into PTM Build Pipeline

**Files:**
- Modify: `hvantk/ptm/pipeline.py`
- Modify: `hvantk/commands/ptm_cli.py`

- [ ] **Step 1: Add `cptac_tsv` to PTMBuildConfig**

In `hvantk/ptm/pipeline.py`, add after `peptideatlas_tsv` field (line 58):

```python
    cptac_tsv: Optional[str] = None
```

Add validation after the PeptideAtlas check (line 75):

```python
        if self.cptac_tsv and not os.path.exists(self.cptac_tsv):
            errors.append(f"CPTAC TSV file not found: {self.cptac_tsv}")
```

- [ ] **Step 2: Add CPTAC concatenation in ptm_build_pipeline()**

In `hvantk/ptm/pipeline.py`, after the PeptideAtlas concatenation block (around line 347, after the `logger.info("Combined:...")` line), add:

```python
    # Step 4c: Map CPTAC sites (if provided)
    if config.cptac_tsv:
        logger.info("Mapping CPTAC phospho sites...")
        cptac_mapped_path = os.path.join(config.output_dir, "cptac_sites_mapped.tsv")
        cptac_result = map_ptm_sites(config.cptac_tsv, gtf_data, cptac_mapped_path)

        # Concatenate with existing mapped TSV
        combined_path = os.path.join(config.output_dir, "ptm_sites_combined.tsv")
        with open(combined_path, "w", newline="") as fout:
            writer = csv.DictWriter(
                fout, fieldnames=PTM_OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for src_path in [mapped_path, cptac_mapped_path]:
                with open(src_path) as fin:
                    reader = csv.DictReader(fin, delimiter="\t")
                    for row in reader:
                        writer.writerow(row)

        mapped_path = combined_path
        result.n_total += cptac_result.n_total
        result.n_mapped += cptac_result.n_mapped
        result.n_failed += cptac_result.n_failed
        for method, count in cptac_result.resolution_counts.items():
            result.resolution_counts[method] = (
                result.resolution_counts.get(method, 0) + count
            )

        logger.info(
            f"Combined: {result.n_mapped} total mapped sites "
            f"(including CPTAC)"
        )
```

- [ ] **Step 3: Add `--cptac-tsv` option to ptm build CLI**

In `hvantk/commands/ptm_cli.py`, add after the `--peptideatlas-tsv` option block (around line 85):

```python
@click.option(
    "--cptac-tsv",
    type=click.Path(exists=True),
    default=None,
    help="Path to CPTAC phospho intermediate TSV (from 'hvantk download cptac-phospho')",
)
```

Update the function signature to include `cptac_tsv`:

```python
def ptm_build(
    ctx, output_dir, output_ht, gtf_path, ptm_tsv, peptideatlas_tsv, cptac_tsv, flanking_codons, overwrite
):
```

Add it to the PTMBuildConfig constructor:

```python
        config = PTMBuildConfig(
            output_dir=output_dir,
            output_ht=output_ht,
            gtf_path=gtf_path,
            ptm_tsv=ptm_tsv,
            peptideatlas_tsv=peptideatlas_tsv,
            cptac_tsv=cptac_tsv,
            flanking_codons=flanking_codons,
            overwrite=overwrite,
        )
```

- [ ] **Step 4: Run tests**

Run: `pytest hvantk/tests/test_ptm.py hvantk/tests/test_cptac_phospho.py -v`
Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add hvantk/ptm/pipeline.py hvantk/commands/ptm_cli.py
git commit -m "feat: integrate CPTAC phospho source into PTM build pipeline"
```

---

### Task 5: Add MatrixTable Builder and CLI

**Files:**
- Modify: `hvantk/tables/matrix_builders.py`
- Modify: `hvantk/tables/registry.py`
- Modify: `hvantk/commands/make_matrix_cli.py`

- [ ] **Step 1: Add `build_cptac_phospho_mt()` to matrix_builders.py**

In `hvantk/tables/matrix_builders.py`, after the `build_cptac_mt()` function (around line 204), add:

```python
def build_cptac_phospho_mt(
    expression_path: str,
    metadata_path: str,
    output_mt: Optional[str] = None,
    site_id_col: str = "SiteID",
    sample_id_col: str = "SampleID",
    categorical_cols: Optional[List[str]] = None,
    numeric_cols: Optional[List[str]] = None,
    overwrite: bool = False,
) -> "hl.MatrixTable":
    """Build a CPTAC phospho MatrixTable from a sites-by-samples matrix.

    Parameters
    ----------
    expression_path : str
        Path to matrix CSV (sites × samples, from cptac_phospho_datasets).
    metadata_path : str
        Path to metadata CSV (sample clinical info).
    output_mt : str, optional
        Path to checkpoint the MatrixTable.
    site_id_col : str
        Column name for site identifiers (default: "SiteID").
    sample_id_col : str
        Column name for sample identifiers in metadata (default: "SampleID").
    categorical_cols : list of str, optional
        Metadata columns to cast as string.
    numeric_cols : list of str, optional
        Metadata columns to cast as float.
    overwrite : bool
        If True, overwrite existing output.

    Returns
    -------
    hl.MatrixTable
    """
    import hail as hl

    logger.info("Building CPTAC phospho MatrixTable")

    # Read expression matrix (sites × samples)
    expr_df = pd.read_csv(expression_path, index_col=0)
    logger.info("Expression matrix: %d sites x %d samples", *expr_df.shape)

    # Melt to coordinate format: SiteID, SampleID, Intensity
    expr_long = expr_df.stack(dropna=False).reset_index()
    expr_long.columns = [site_id_col, sample_id_col, "Intensity"]
    expr_long[site_id_col] = expr_long[site_id_col].astype(str)
    expr_long[sample_id_col] = expr_long[sample_id_col].astype(str)

    # Parse site ID into gene_symbol, amino_acid, residue_pos
    def _parse_site_id(sid):
        parts = sid.rsplit("_", 1)
        if len(parts) == 2:
            gene = parts[0]
            site = parts[1]
            if site and site[0] in "STY" and site[1:].isdigit():
                return gene, site[0], int(site[1:])
        return sid, "", 0

    site_info = {sid: _parse_site_id(sid) for sid in expr_long[site_id_col].unique()}
    expr_long["gene_symbol"] = expr_long[site_id_col].map(lambda s: site_info[s][0])
    expr_long["amino_acid"] = expr_long[site_id_col].map(lambda s: site_info[s][1])
    expr_long["residue_pos"] = expr_long[site_id_col].map(lambda s: site_info[s][2])

    # Create Hail Table from long-format expression
    ht_expr = hl.Table.from_pandas(expr_long)
    ht_expr = ht_expr.key_by(site_id_col, sample_id_col)

    # Convert to MatrixTable
    mt = ht_expr.to_matrix_table(
        row_key=[site_id_col],
        col_key=[sample_id_col],
        row_fields=["gene_symbol", "amino_acid", "residue_pos"],
    )

    # Read and join metadata
    meta_df = pd.read_csv(metadata_path, index_col=0)
    meta_df.index = meta_df.index.astype(str)
    meta_df.index.name = sample_id_col

    if categorical_cols:
        for c in categorical_cols:
            if c in meta_df.columns:
                meta_df[c] = meta_df[c].astype(str)
    if numeric_cols:
        for c in numeric_cols:
            if c in meta_df.columns:
                meta_df[c] = pd.to_numeric(meta_df[c], errors="coerce")

    ht_meta = hl.Table.from_pandas(meta_df.reset_index())
    ht_meta = ht_meta.key_by(sample_id_col)
    mt = mt.annotate_cols(**ht_meta[mt.col_key])

    mt = annotate_column_summary(mt)
    mt = mt.annotate_globals(
        hvantk_metadata=build_matrix_metadata("CPTAC-phospho", expression_path, mt)
    )

    if output_mt:
        logger.info("Checkpointing CPTAC phospho MatrixTable to %s", output_mt)
        mt = mt.checkpoint(output_mt, overwrite=overwrite)

    return mt
```

Also ensure the necessary imports exist at the top of the file. Add if needed:

```python
import pandas as pd
```

And ensure `annotate_column_summary` and `build_matrix_metadata` are imported from the appropriate module.

- [ ] **Step 2: Register in MATRIX_BUILDERS**

In `hvantk/tables/registry.py`, add after the existing `"cptac"` entry in MATRIX_BUILDERS (around line 311):

```python
    "cptac-phospho": create_matrix_adapter(
        "hvantk.tables.matrix_builders",
        "build_cptac_phospho_mt",
        required_inputs=["expression", "metadata"],
    ),
```

- [ ] **Step 3: Add mkmatrix cptac-phospho CLI command**

In `hvantk/commands/make_matrix_cli.py`, after the existing `mkmatrix_cptac` command (around line 236), add:

```python
@mkmatrix_group.command("cptac-phospho")
@click.option(
    "-e",
    "--expression",
    required=True,
    type=click.Path(exists=True),
    help="Path to CPTAC phospho matrix CSV (sites × samples)",
)
@click.option(
    "-m",
    "--metadata",
    required=True,
    type=click.Path(exists=True),
    help="Path to CPTAC phospho metadata CSV",
)
@click.option("-o", "--output-mt", "output_mt", required=True, type=click.Path())
@click.option("--sid", "--sample-id-col", "sample_id_col", default="SampleID", show_default=True)
@click.option(
    "-c",
    "--categorical-cols",
    default=None,
    help="Comma-separated categorical metadata columns",
)
@click.option(
    "-u",
    "--numeric-cols",
    default=None,
    help="Comma-separated numeric metadata columns",
)
@click.option("-w", "--overwrite", is_flag=True)
def mkmatrix_cptac_phospho(
    expression,
    metadata,
    output_mt,
    sample_id_col,
    categorical_cols,
    numeric_cols,
    overwrite,
):
    """Build a CPTAC phospho MatrixTable from site intensities and metadata."""
    logger.info("Building CPTAC phospho MatrixTable")

    def _split(val):
        if not val:
            return None
        parts = [x.strip() for x in val.split(",") if x.strip()]
        return parts or None

    from hvantk.tables.matrix_builders import build_cptac_phospho_mt

    mt = build_cptac_phospho_mt(
        expression_path=expression,
        metadata_path=metadata,
        output_mt=output_mt,
        sample_id_col=sample_id_col,
        categorical_cols=_split(categorical_cols),
        numeric_cols=_split(numeric_cols),
        overwrite=overwrite,
    )
    click.echo(f"MatrixTable created at {output_mt}")
    mt.describe()
```

- [ ] **Step 4: Verify CLI registration**

Run: `python -c "from hvantk.commands.make_matrix_cli import mkmatrix_group; from click.testing import CliRunner; r = CliRunner().invoke(mkmatrix_group, ['cptac-phospho', '--help']); print(r.output); assert '--expression' in r.output"`

- [ ] **Step 5: Commit**

```bash
git add hvantk/tables/matrix_builders.py hvantk/tables/registry.py hvantk/commands/make_matrix_cli.py
git commit -m "feat: add cptac-phospho MatrixTable builder and CLI command"
```

---

### Task 6: Final Verification

- [ ] **Step 1: Run full test suite**

Run: `pytest hvantk/tests/test_ptm.py hvantk/tests/test_cptac_phospho.py hvantk/tests/test_peptideatlas_phospho.py -v --tb=short`
Expected: ALL PASS

- [ ] **Step 2: Verify all CLI commands**

```python
python -c "
from click.testing import CliRunner
from hvantk.commands.download_cli import download_group
from hvantk.commands.ptm_cli import ptm_group
from hvantk.commands.make_matrix_cli import mkmatrix_group

runner = CliRunner()

# Download help
r = runner.invoke(download_group, ['--help'])
assert 'cptac-phospho' in r.output

# List cancers
r = runner.invoke(download_group, ['cptac-phospho', '--list-cancers'])
assert 'brca' in r.output

# PTM build help
r = runner.invoke(ptm_group, ['build', '--help'])
assert '--cptac-tsv' in r.output

# mkmatrix cptac-phospho help
r = runner.invoke(mkmatrix_group, ['cptac-phospho', '--help'])
assert '--expression' in r.output

print('All CLI checks passed!')
"
```

- [ ] **Step 3: Commit any final fixes**

If any tests failed or CLI output was wrong, fix and commit.
