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
    """Mock phosphoproteomics DataFrame mimicking cptac package output."""
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
            [1.5, 2.0, 0.5, np.nan],
            [np.nan, 1.8, -0.3, 3.2],
            [0.7, np.nan, 1.1, 2.5],
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

    # MAPK1 T185_Y187 -> two separate sites (T185 and Y187)
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
    assert set(df.columns) == {"Sample_01", "Sample_02", "Sample_03"}
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
