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


def test_aggregate_by_disease_category(clingen_table_path):
    streamer = ClinGenStreamer(clingen_table_path)
    categories = {"cancer": ["cancer"]}
    result = streamer.aggregate_by_disease_category(categories)
    assert "cancer" in result
    assert "BRCA1" in result["cancer"]
