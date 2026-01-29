"""
Tests for gene set data structures and I/O.
"""

import json
import pytest
from pathlib import Path

from hvantk.enrichex.gene_sets import (
    GeneSet,
    GeneSetCollection,
    load_gene_sets_from_dict,
    load_marker_genes,
)


class TestGeneSet:
    """Tests for GeneSet class."""

    def test_gene_set_creation(self):
        """Test basic GeneSet creation."""
        gs = GeneSet(
            name="test_set",
            genes={"BRCA1", "TP53", "EGFR"},
            source="test",
            metadata={"method": "manual"},
        )

        assert gs.name == "test_set"
        assert gs.n_genes == 3
        assert "BRCA1" in gs.genes
        assert gs.source == "test"
        assert gs.metadata["method"] == "manual"

    def test_gene_set_to_dict(self):
        """Test GeneSet serialization to dictionary."""
        gs = GeneSet(
            name="test_set",
            genes={"BRCA1", "TP53"},
            source="test",
        )

        d = gs.to_dict()

        assert d["name"] == "test_set"
        assert d["n_genes"] == 2
        assert isinstance(d["genes"], list)
        assert set(d["genes"]) == {"BRCA1", "TP53"}
        # Genes should be sorted for reproducibility
        assert d["genes"] == sorted(d["genes"])

    def test_gene_set_from_dict(self):
        """Test GeneSet deserialization from dictionary."""
        data = {
            "name": "test_set",
            "genes": ["BRCA1", "TP53"],
            "source": "test",
            "metadata": {"key": "value"},
        }

        gs = GeneSet.from_dict(data)

        assert gs.name == "test_set"
        assert gs.n_genes == 2
        assert gs.genes == {"BRCA1", "TP53"}
        assert gs.source == "test"
        assert gs.metadata["key"] == "value"

    def test_gene_set_round_trip(self):
        """Test GeneSet serialization round-trip."""
        gs1 = GeneSet(
            name="test",
            genes={"A", "B", "C"},
            source="manual",
            metadata={"x": 1},
        )

        # Round-trip through dict
        d = gs1.to_dict()
        gs2 = GeneSet.from_dict(d)

        assert gs1.name == gs2.name
        assert gs1.genes == gs2.genes
        assert gs1.source == gs2.source
        assert gs1.metadata == gs2.metadata


class TestGeneSetCollection:
    """Tests for GeneSetCollection class."""

    def test_gene_set_collection_creation(self):
        """Test basic GeneSetCollection creation."""
        gs1 = GeneSet("set1", {"A", "B", "C"})
        gs2 = GeneSet("set2", {"C", "D", "E"})

        collection = GeneSetCollection(
            gene_sets={"set1": gs1, "set2": gs2},
            background_genes={"A", "B", "C", "D", "E", "F"},
            source_description="test",
        )

        assert len(collection) == 2
        assert "set1" in collection.gene_sets
        assert "set2" in collection.gene_sets
        assert len(collection.background_genes) == 6

    def test_gene_set_collection_iteration(self):
        """Test iterating over GeneSetCollection."""
        gs1 = GeneSet("set1", {"A", "B"})
        gs2 = GeneSet("set2", {"C", "D"})

        collection = GeneSetCollection(
            gene_sets={"set1": gs1, "set2": gs2},
            background_genes={"A", "B", "C", "D"},
        )

        gene_sets = list(collection)
        assert len(gene_sets) == 2
        assert gs1 in gene_sets
        assert gs2 in gene_sets

    def test_gene_set_collection_get(self):
        """Test getting gene sets by name."""
        gs = GeneSet("test", {"A", "B"})
        collection = GeneSetCollection(
            gene_sets={"test": gs}, background_genes={"A", "B"}
        )

        assert collection.get("test") == gs
        assert collection.get("nonexistent") is None

    def test_gene_set_collection_names(self):
        """Test getting gene set names."""
        gs1 = GeneSet("zebra", {"A"})
        gs2 = GeneSet("alpha", {"B"})

        collection = GeneSetCollection(
            gene_sets={"zebra": gs1, "alpha": gs2},
            background_genes={"A", "B"},
        )

        names = collection.names()
        # Should be sorted
        assert names == ["alpha", "zebra"]

    def test_gene_set_collection_to_dict(self):
        """Test GeneSetCollection serialization."""
        gs1 = GeneSet("set1", {"A", "B"})
        collection = GeneSetCollection(
            gene_sets={"set1": gs1},
            background_genes={"A", "B", "C"},
            source_description="test",
        )

        d = collection.to_dict()

        assert "gene_sets" in d
        assert "background_genes" in d
        assert d["n_gene_sets"] == 1
        assert d["n_background"] == 3
        assert d["source_description"] == "test"

    def test_gene_set_collection_from_dict(self):
        """Test GeneSetCollection deserialization."""
        data = {
            "gene_sets": {
                "set1": {
                    "name": "set1",
                    "genes": ["A", "B"],
                    "source": "test",
                    "metadata": {},
                }
            },
            "background_genes": ["A", "B", "C"],
            "source_description": "test",
            "metadata": {},
        }

        collection = GeneSetCollection.from_dict(data)

        assert len(collection) == 1
        assert collection.get("set1") is not None
        assert len(collection.background_genes) == 3

    def test_gene_set_collection_save_load(self, tmp_path):
        """Test saving and loading GeneSetCollection."""
        gs1 = GeneSet("set1", {"A", "B", "C"})
        gs2 = GeneSet("set2", {"D", "E", "F"})

        collection = GeneSetCollection(
            gene_sets={"set1": gs1, "set2": gs2},
            background_genes={"A", "B", "C", "D", "E", "F"},
            source_description="test collection",
        )

        # Save
        path = tmp_path / "gene_sets.json"
        collection.save(path)

        # Check file exists and is valid JSON
        assert path.exists()
        with open(path) as f:
            data = json.load(f)
        assert "gene_sets" in data

        # Load
        loaded = GeneSetCollection.load(path)

        assert len(loaded) == len(collection)
        assert loaded.get("set1").genes == gs1.genes
        assert loaded.get("set2").genes == gs2.genes
        assert loaded.background_genes == collection.background_genes


class TestGeneSetLoadingFunctions:
    """Tests for gene set loading utility functions."""

    def test_load_gene_sets_from_dict(self):
        """Test creating GeneSetCollection from simple dict."""
        gene_sets_dict = {
            "set1": ["A", "B", "C"],
            "set2": ["D", "E", "F"],
        }

        collection = load_gene_sets_from_dict(gene_sets_dict, source="test")

        assert len(collection) == 2
        assert collection.get("set1").n_genes == 3
        assert collection.get("set2").n_genes == 3
        # Background should be union of all genes
        assert len(collection.background_genes) == 6

    def test_load_gene_sets_from_dict_with_background(self):
        """Test creating GeneSetCollection with custom background."""
        gene_sets_dict = {
            "set1": ["A", "B"],
        }
        background = {"A", "B", "C", "D", "E"}

        collection = load_gene_sets_from_dict(
            gene_sets_dict, background_genes=background
        )

        assert len(collection.background_genes) == 5
        assert "E" in collection.background_genes

    def test_load_marker_genes(self):
        """Test loading marker genes from TSV file."""
        # Use test data
        test_file = Path(__file__).parent / "testdata" / "marker_genes.tsv"

        if not test_file.exists():
            pytest.skip("Test data file not found")

        collection = load_marker_genes(
            test_file, cluster_column="cluster", gene_column="gene"
        )

        # Should have 4 clusters: T_cell, B_cell, Macrophage, NK_cell
        assert len(collection) == 4
        assert collection.get("T_cell") is not None
        assert collection.get("B_cell") is not None

        # T_cell should have specific genes
        t_cell_genes = collection.get("T_cell").genes
        assert "CD3D" in t_cell_genes
        assert "CD8A" in t_cell_genes

    def test_load_marker_genes_missing_column(self, tmp_path):
        """Test error handling for missing column."""
        # Create invalid TSV
        invalid_file = tmp_path / "invalid.tsv"
        invalid_file.write_text("wrong_column\tgene\nvalue\tGENE1\n")

        with pytest.raises(ValueError, match="Column 'cluster' not found"):
            load_marker_genes(invalid_file, cluster_column="cluster")


class TestGeneSetIntegration:
    """Integration tests for gene set workflows."""

    def test_complete_workflow(self, tmp_path):
        """Test complete workflow: create, save, load, use."""
        # Create gene sets
        gene_sets = {
            "cancer_genes": ["BRCA1", "TP53", "EGFR"],
            "cell_cycle": ["CDK1", "CDK2", "CCNA1"],
        }

        # Create collection
        collection = load_gene_sets_from_dict(gene_sets, source="test")

        # Save
        path = tmp_path / "gene_sets.json"
        collection.save(path)

        # Load
        loaded = GeneSetCollection.load(path)

        # Verify
        assert len(loaded) == 2
        assert loaded.get("cancer_genes").n_genes == 3
        assert "BRCA1" in loaded.get("cancer_genes").genes
