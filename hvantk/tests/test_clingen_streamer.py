"""
Hail integration tests for ClinGenStreamer.
"""

from pathlib import Path

import pytest

from hvantk.data.clingen_streamer import ClinGenStreamer
from hvantk.tables.table_builders import create_clingen_gene_disease_tb

pytestmark = [pytest.mark.hail, pytest.mark.slow]

TEST_DIR = Path(__file__).parent / "testdata"


@pytest.fixture
def clingen_table_path(tmp_path):
    input_path = TEST_DIR / "raw/clingen/clingen_test_sample.csv"
    output_path = tmp_path / "clingen_test.ht"
    create_clingen_gene_disease_tb(
        input_path=str(input_path),
        output_path=str(output_path),
        overwrite=True,
    )
    return str(output_path)


def test_get_genes_by_classification_definitive(clingen_table_path):
    streamer = ClinGenStreamer(clingen_table_path)
    genes = streamer.get_genes_by_classification("Definitive")
    assert "BRCA1" in genes
    assert "APOB" not in genes


def test_get_genes_by_disease_cancer(clingen_table_path):
    streamer = ClinGenStreamer(clingen_table_path)
    genes = streamer.get_genes_by_disease("cancer", match_mode="contains")
    assert "BRCA1" in genes
    assert "BRCA2" in genes


def test_get_genes_by_mondo_id(clingen_table_path):
    streamer = ClinGenStreamer(clingen_table_path)
    genes = streamer.get_genes_by_mondo_id("MONDO:0005144")
    assert genes == {"BRCA1"}


def test_compute_stats(clingen_table_path):
    streamer = ClinGenStreamer(clingen_table_path)
    stats = streamer.compute_stats()
    assert stats["total_associations"] == 11
    assert stats["unique_genes"] == 10
    assert "Definitive" in stats["classification_counts"]


def test_get_geneset_per_gcep(clingen_table_path):
    streamer = ClinGenStreamer(clingen_table_path)
    result = streamer.get_geneset_per_gcep()
    # Test data has 4 distinct GCEPs
    assert len(result) >= 3
    # "Hereditary Cancer GCEP" should be shortened to "Hereditary Cancer"
    assert "Hereditary Cancer" in result
    assert "TP53" in result["Hereditary Cancer"]
    assert "PTEN" in result["Hereditary Cancer"]
    assert "CDH1" in result["Hereditary Cancer"]
    # Breast/Ovarian panel
    hereditary_bop = "Hereditary Breast, Ovarian and Pancreatic Cancer"
    assert hereditary_bop in result
    assert "BRCA1" in result[hereditary_bop]
    assert "BRCA2" in result[hereditary_bop]


def test_get_geneset_per_gcep_min_genes(clingen_table_path):
    streamer = ClinGenStreamer(clingen_table_path)
    # Require at least 3 genes per GCEP — should filter out small panels
    result = streamer.get_geneset_per_gcep(min_genes=3)
    for genes in result.values():
        assert len(genes) >= 3


def test_get_geneset_per_gcep_min_classification(clingen_table_path):
    streamer = ClinGenStreamer(clingen_table_path)
    result_all = streamer.get_geneset_per_gcep()
    result_definitive = streamer.get_geneset_per_gcep(min_classification="Definitive")
    # Filtering to Definitive should yield fewer or equal genes
    for gcep in result_definitive:
        if gcep in result_all:
            assert result_definitive[gcep] <= result_all[gcep]


def test_get_geneset_per_gcep_shorten_names_false(clingen_table_path):
    streamer = ClinGenStreamer(clingen_table_path)
    result = streamer.get_geneset_per_gcep(shorten_names=False)
    # Should contain the full GCEP name
    assert any("GCEP" in name for name in result)


def test_aggregate_by_disease_category(clingen_table_path):
    streamer = ClinGenStreamer(clingen_table_path)
    categories = {"cancer": ["cancer"]}
    result = streamer.aggregate_by_disease_category(categories)
    assert "cancer" in result
    assert "BRCA1" in result["cancer"]
