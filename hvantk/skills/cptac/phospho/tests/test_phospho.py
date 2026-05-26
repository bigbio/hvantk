"""Unit tests for CPTAC phospho parser helpers and lifecycle entry points.

These tests exercise the shared helper functions (``parse_phospho_site``,
``extract_phospho_sites``, ``write_intermediate_tsv``, ``parse_raw_dir``)
without requiring the upstream ``cptac`` Python package. Builder coverage
lives in :mod:`hvantk.tests.test_expression_builders_anndata`.
"""

import csv
import os

import pandas as pd
import pytest


def test_parse_phospho_site_single():
    from hvantk.skills.cptac.shared.datasets import parse_phospho_site

    assert parse_phospho_site("S65") == [("S", 65)]
    assert parse_phospho_site("T185") == [("T", 185)]
    assert parse_phospho_site("Y243") == [("Y", 243)]


def test_parse_phospho_site_multi():
    from hvantk.skills.cptac.shared.datasets import parse_phospho_site

    assert parse_phospho_site("T185_Y187") == [("T", 185), ("Y", 187)]


def test_parse_phospho_site_unknown_residue():
    from hvantk.skills.cptac.shared.datasets import parse_phospho_site

    # K (lysine) is not a phospho residue, should be filtered out
    assert parse_phospho_site("K12") == []
    assert parse_phospho_site("") == []


def test_extract_phospho_sites_aggregates_duplicates():
    from hvantk.skills.cptac.shared.datasets import extract_phospho_sites

    columns = pd.MultiIndex.from_tuples(
        [
            ("TP53", "S15", "peptideA", "P04637"),
            ("TP53", "S15", "peptideB", "P04637"),
            ("EGFR", "Y1068", "peptideC", "P00533"),
        ],
        names=["Gene", "Site", "Peptide", "Database_ID"],
    )
    df = pd.DataFrame(
        {
            columns[0]: [1.0, 2.0, 3.0],
            columns[1]: [4.0, 5.0, 6.0],
            columns[2]: [7.0, 8.0, 9.0],
        }
    )
    sites = extract_phospho_sites(df, cancer_type="brca", tissue_type="tumor")

    assert len(sites) == 2  # TP53_S15 (dedup), EGFR_Y1068
    tp53 = [s for s in sites if s["gene_symbol"] == "TP53"][0]
    egfr = [s for s in sites if s["gene_symbol"] == "EGFR"][0]
    assert tp53["amino_acid"] == "S"
    assert tp53["position"] == 15
    assert tp53["cancer_type"] == "brca"
    assert tp53["tissue_type"] == "tumor"
    assert tp53["source_db"] == "CPTAC"
    assert tp53["n_observations"] == 6  # 3 + 3
    assert egfr["amino_acid"] == "Y"
    assert egfr["position"] == 1068


def test_write_intermediate_tsv_round_trip(tmp_path):
    from hvantk.skills.cptac.shared.datasets import (
        _TSV_COLUMNS,
        write_intermediate_tsv,
    )

    sites = [
        {
            "accession": "P04637",
            "gene_symbol": "TP53",
            "position": 15,
            "description": "Phosphoserine",
            "amino_acid": "S",
            "ensembl_xrefs": "",
            "sequence_length": "",
            "n_observations": 6,
            "source_db": "CPTAC",
            "evidence_type": "mass_spectrometry",
            "tissue_type": "tumor",
            "cancer_type": "brca",
            "mean_intensity": 1.5,
        }
    ]
    out_path = tmp_path / "out.tsv"
    write_intermediate_tsv(sites, str(out_path))

    with open(out_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["gene_symbol"] == "TP53"
    assert rows[0]["position"] == "15"
    assert set(_TSV_COLUMNS).issubset(rows[0].keys())


def test_parse_raw_dir_consolidates_per_cancer_tsvs(tmp_path):
    from hvantk.skills.cptac.shared.datasets import (
        _TSV_COLUMNS,
        parse_raw_dir,
        write_intermediate_tsv,
    )

    brca_sites = [
        {
            "accession": "P04637",
            "gene_symbol": "TP53",
            "position": 15,
            "description": "Phosphoserine",
            "amino_acid": "S",
            "n_observations": 6,
            "source_db": "CPTAC",
            "evidence_type": "mass_spectrometry",
            "tissue_type": "tumor",
            "cancer_type": "brca",
            "mean_intensity": 1.5,
        }
    ]
    luad_sites = [
        {
            "accession": "P00533",
            "gene_symbol": "EGFR",
            "position": 1068,
            "description": "Phosphotyrosine",
            "amino_acid": "Y",
            "n_observations": 3,
            "source_db": "CPTAC",
            "evidence_type": "mass_spectrometry",
            "tissue_type": "tumor",
            "cancer_type": "luad",
            "mean_intensity": 2.5,
        }
    ]
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    write_intermediate_tsv(brca_sites, str(raw_dir / "cptac-phospho-brca.tsv"))
    write_intermediate_tsv(luad_sites, str(raw_dir / "cptac-phospho-luad.tsv"))
    # Tumor-only files should be ignored
    write_intermediate_tsv(brca_sites, str(raw_dir / "cptac-phospho-brca-tumor.tsv"))

    out_path = tmp_path / "merged.tsv"
    result = parse_raw_dir(str(raw_dir), str(out_path))
    assert result == str(out_path)

    with open(out_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    genes = sorted(r["gene_symbol"] for r in rows)
    assert genes == ["EGFR", "TP53"]


def test_parse_raw_dir_filters_by_cancer_type(tmp_path):
    from hvantk.skills.cptac.shared.datasets import (
        parse_raw_dir,
        write_intermediate_tsv,
    )

    brca_sites = [
        {
            "accession": "P04637",
            "gene_symbol": "TP53",
            "position": 15,
            "amino_acid": "S",
            "n_observations": 6,
            "tissue_type": "tumor",
            "cancer_type": "brca",
            "mean_intensity": 1.5,
        }
    ]
    luad_sites = [
        {
            "accession": "P00533",
            "gene_symbol": "EGFR",
            "position": 1068,
            "amino_acid": "Y",
            "n_observations": 3,
            "tissue_type": "tumor",
            "cancer_type": "luad",
            "mean_intensity": 2.5,
        }
    ]
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    write_intermediate_tsv(brca_sites, str(raw_dir / "cptac-phospho-brca.tsv"))
    write_intermediate_tsv(luad_sites, str(raw_dir / "cptac-phospho-luad.tsv"))

    out_path = tmp_path / "brca.tsv"
    parse_raw_dir(str(raw_dir), str(out_path), cancer_type="brca")

    with open(out_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["gene_symbol"] == "TP53"


def test_parse_raw_dir_missing_input_raises(tmp_path):
    from hvantk.skills.cptac.shared.datasets import parse_raw_dir

    with pytest.raises(FileNotFoundError):
        parse_raw_dir(str(tmp_path), str(tmp_path / "out.tsv"))


def test_dataset_class_rejects_unknown_cancer_type():
    from hvantk.skills.cptac.shared.datasets import CPTACPhosphoDataset

    with pytest.raises(ValueError):
        CPTACPhosphoDataset(cancer_type="not-a-cancer")


def test_dataset_class_metadata_includes_class_name():
    from hvantk.skills.cptac.shared.datasets import CPTACPhosphoDataset

    ds = CPTACPhosphoDataset(cancer_type="brca")
    meta = ds.get_metadata()
    assert meta["cancer_type"] == "brca"
    assert meta["class_name"] == "Brca"
