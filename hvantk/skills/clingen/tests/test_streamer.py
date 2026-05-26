"""
Hail integration tests for ClinGenStreamer.
"""

from pathlib import Path

import pytest

from hvantk.skills.clingen.streamer import ClinGenStreamer
from hvantk.skills.clingen.builder import create_clingen_gene_disease_tb

pytestmark = [pytest.mark.hail, pytest.mark.slow]

# Fixture lives under the clingen plugin folder after the Phase 2 migration.
TEST_FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "skills"
    / "clingen"
    / "tests"
    / "testdata"
    / "raw"
    / "clingen"
    / "clingen_test_sample.csv"
)


@pytest.fixture
def clingen_table_path(tmp_path):
    input_path = TEST_FIXTURE
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


def test_get_genes_by_classification_with_gene_mapper(clingen_table_path, tmp_path):
    """Test GeneMapper integration translates gene symbols to HGNC IDs."""
    from unittest.mock import Mock

    streamer = ClinGenStreamer(clingen_table_path)

    # Get baseline symbols
    symbols = streamer.get_genes_by_classification("Definitive")
    assert len(symbols) > 0

    # Mock GeneMapper to verify it's called correctly
    mock_mapper = Mock()
    mock_mapper.map_ids.return_value = {s: f"HGNC:{i}" for i, s in enumerate(symbols)}

    result = streamer.get_genes_by_classification(
        "Definitive",
        gene_mapper=mock_mapper,
        output_id_type="hgnc_id",
    )

    mock_mapper.map_ids.assert_called_once()
    call_args = mock_mapper.map_ids.call_args
    assert call_args.kwargs["source_type"] == "gene_symbol"
    assert call_args.kwargs["target_type"] == "hgnc_id"
    # Results should be the mapped HGNC IDs
    assert all(v.startswith("HGNC:") for v in result)


def test_to_gene_set_with_gene_mapper(clingen_table_path):
    """Test to_gene_set with GeneMapper translates IDs."""
    from unittest.mock import Mock

    streamer = ClinGenStreamer(clingen_table_path)

    mock_mapper = Mock()
    mock_mapper.map_ids.return_value = {"BRCA1": "ENSG00000012048"}

    result = streamer.to_gene_set(
        min_classification="Definitive",
        gene_mapper=mock_mapper,
        output_id_type="ensembl_gene_id",
    )

    mock_mapper.map_ids.assert_called_once()
    assert "ENSG00000012048" in result


def test_gene_mapper_none_does_not_translate(clingen_table_path):
    """Test that without gene_mapper, symbols are returned as-is."""
    streamer = ClinGenStreamer(clingen_table_path)
    result = streamer.get_genes_by_classification("Definitive")
    # Should be gene symbols, not HGNC IDs or Ensembl IDs
    assert all(not v.startswith("HGNC:") for v in result)
    assert all(not v.startswith("ENSG") for v in result)


def test_get_genes_by_classification_output_id_requires_gene_mapper(
    clingen_table_path,
):
    """output_id_type without gene_mapper raises ValueError."""
    streamer = ClinGenStreamer(clingen_table_path)
    with pytest.raises(ValueError, match="output_id_type was provided"):
        streamer.get_genes_by_classification(
            "Definitive",
            output_id_type="hgnc_id",
        )


def test_to_gene_set_output_id_requires_gene_mapper(clingen_table_path):
    """to_gene_set requires gene_mapper when output_id_type is used."""
    streamer = ClinGenStreamer(clingen_table_path)
    with pytest.raises(ValueError, match="output_id_type was provided"):
        streamer.to_gene_set(
            min_classification="Definitive",
            output_id_type="hgnc_id",
        )
