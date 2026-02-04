"""
Hail integration tests for HGNC gene nomenclature table builder and GeneMapper.

These tests require Hail and use the test fixture TSV file.
"""

import pytest
import shutil
from pathlib import Path

from hvantk.tables.table_builders import create_hgnc_gene_tb
from hvantk.data.gene_mapper import GeneMapper

# Mark as Hail-dependent and slow
pytestmark = [pytest.mark.hail, pytest.mark.slow]

# Test data directory
TEST_DIR = Path(__file__).parent / "testdata"
# Temporary directory for testing
TMP_DIR = Path(__file__).parent / "tmp"


@pytest.fixture(autouse=True)
def setup_teardown():
    """Create and clean up the temporary directory for each test."""
    TMP_DIR.mkdir(exist_ok=True, parents=True)
    yield
    if TMP_DIR.exists():
        shutil.rmtree(TMP_DIR)


class TestCreateHgncGeneTb:
    """Tests for the create_hgnc_gene_tb function."""

    def test_create_hgnc_gene_tb_default(self):
        """Test building HGNC table with default settings."""
        input_path = TEST_DIR / "raw/hgnc/hgnc_test_sample.tsv"
        output_path = TMP_DIR / "hgnc.ht"

        tb = create_hgnc_gene_tb(
            input_path=str(input_path),
            output_path=str(output_path),
            overwrite=True,
        )

        # Check row count (5 approved genes, 1 withdrawn excluded)
        assert tb.count() == 5

        # Check key field
        key_fields = list(tb.key.dtype)
        assert key_fields == ["hgnc_id"]

        # Check that _SUCCESS file exists
        success_file = Path(output_path) / "_SUCCESS"
        assert success_file.exists()

    def test_create_hgnc_gene_tb_include_withdrawn(self):
        """Test building HGNC table including withdrawn genes."""
        input_path = TEST_DIR / "raw/hgnc/hgnc_test_sample.tsv"
        output_path = TMP_DIR / "hgnc_with_withdrawn.ht"

        tb = create_hgnc_gene_tb(
            input_path=str(input_path),
            output_path=str(output_path),
            include_withdrawn=True,
            overwrite=True,
        )

        # Check row count (all 6 genes including withdrawn)
        assert tb.count() == 6

    def test_create_hgnc_gene_tb_field_mapping(self):
        """Test that fields are correctly renamed."""
        input_path = TEST_DIR / "raw/hgnc/hgnc_test_sample.tsv"
        output_path = TMP_DIR / "hgnc_fields.ht"

        tb = create_hgnc_gene_tb(
            input_path=str(input_path),
            output_path=str(output_path),
            overwrite=True,
        )

        # Check renamed fields exist
        row_fields = set(tb.row.dtype)
        assert "hgnc_id" in row_fields
        assert "gene_symbol" in row_fields  # renamed from 'symbol'
        assert "gene_name" in row_fields  # renamed from 'name'
        assert "ensembl_gene_id" in row_fields
        assert "entrez_id" in row_fields
        assert "uniprot_ids" in row_fields

        # Original field names should not exist
        assert "symbol" not in row_fields
        assert "name" not in row_fields

    def test_create_hgnc_gene_tb_pipe_separated_fields(self):
        """Test that pipe-separated fields are parsed into arrays."""
        input_path = TEST_DIR / "raw/hgnc/hgnc_test_sample.tsv"
        output_path = TMP_DIR / "hgnc_arrays.ht"

        tb = create_hgnc_gene_tb(
            input_path=str(input_path),
            output_path=str(output_path),
            overwrite=True,
        )

        # Get BRCA1 row
        brca1 = tb.filter(tb.hgnc_id == "HGNC:1100").collect()[0]

        # Check alias_symbols is an array
        assert isinstance(brca1.alias_symbols, list)
        assert "BRCC1" in brca1.alias_symbols
        assert "FANCS" in brca1.alias_symbols

        # Check uniprot_ids is an array
        assert isinstance(brca1.uniprot_ids, list)
        assert "P38398" in brca1.uniprot_ids

    def test_create_hgnc_gene_tb_empty_arrays(self):
        """Test that empty pipe-separated fields become empty arrays."""
        input_path = TEST_DIR / "raw/hgnc/hgnc_test_sample.tsv"
        output_path = TMP_DIR / "hgnc_empty.ht"

        tb = create_hgnc_gene_tb(
            input_path=str(input_path),
            output_path=str(output_path),
            overwrite=True,
        )

        # GSTT1 has no alias_symbols (empty in test data)
        gstt1 = tb.filter(tb.hgnc_id == "HGNC:4641").collect()[0]
        assert gstt1.alias_symbols == []

    def test_create_hgnc_gene_tb_field_selection(self):
        """Test field selection."""
        input_path = TEST_DIR / "raw/hgnc/hgnc_test_sample.tsv"
        output_path = TMP_DIR / "hgnc_select.ht"

        tb = create_hgnc_gene_tb(
            input_path=str(input_path),
            output_path=str(output_path),
            fields=["gene_symbol", "ensembl_gene_id"],
            overwrite=True,
        )

        row_fields = set(tb.row.dtype)
        # Key field should be present
        assert "hgnc_id" in row_fields
        # Selected fields
        assert "gene_symbol" in row_fields
        assert "ensembl_gene_id" in row_fields
        # Non-selected fields should not be present
        assert "gene_name" not in row_fields
        assert "locus_group" not in row_fields

    def test_create_hgnc_gene_tb_export_tsv(self):
        """Test TSV export functionality."""
        input_path = TEST_DIR / "raw/hgnc/hgnc_test_sample.tsv"
        output_path = TMP_DIR / "hgnc_export.ht"

        tb = create_hgnc_gene_tb(
            input_path=str(input_path),
            output_path=str(output_path),
            overwrite=True,
            export_tsv=True,
        )

        # Check TSV file was created
        tsv_path = Path(str(output_path) + ".tsv.bgz")
        assert tsv_path.exists()


class TestGeneMapper:
    """Tests for the GeneMapper class."""

    @pytest.fixture
    def hgnc_table(self):
        """Create HGNC table for testing."""
        input_path = TEST_DIR / "raw/hgnc/hgnc_test_sample.tsv"
        output_path = TMP_DIR / "hgnc_mapper.ht"

        return create_hgnc_gene_tb(
            input_path=str(input_path),
            output_path=str(output_path),
            overwrite=True,
        )

    @pytest.fixture
    def mapper(self, hgnc_table):
        """Create GeneMapper instance."""
        return GeneMapper(hgnc_table)

    def test_mapper_init(self, hgnc_table):
        """Test GeneMapper initialization."""
        mapper = GeneMapper(hgnc_table)
        assert mapper is not None

    def test_mapper_init_validates_key(self, hgnc_table):
        """Test that mapper validates table key."""
        import hail as hl

        # Re-key table incorrectly
        bad_table = hgnc_table.key_by("gene_symbol")

        with pytest.raises(ValueError, match="keyed by hgnc_id"):
            GeneMapper(bad_table)

    def test_map_to_hgnc_from_symbol(self, mapper):
        """Test mapping gene symbols to HGNC IDs."""
        result = mapper.map_to_hgnc(["BRCA1", "BRCA2", "UNKNOWN"], source_type="gene_symbol")

        assert result["BRCA1"] == "HGNC:1100"
        assert result["BRCA2"] == "HGNC:1101"
        assert result["UNKNOWN"] is None

    def test_map_to_hgnc_from_ensembl(self, mapper):
        """Test mapping Ensembl IDs to HGNC IDs."""
        result = mapper.map_to_hgnc(
            ["ENSG00000012048", "ENSG00000139618"],
            source_type="ensembl_gene_id",
        )

        assert result["ENSG00000012048"] == "HGNC:1100"  # BRCA1
        assert result["ENSG00000139618"] == "HGNC:1101"  # BRCA2

    def test_map_from_hgnc_to_symbol(self, mapper):
        """Test mapping HGNC IDs to gene symbols."""
        result = mapper.map_from_hgnc(
            ["HGNC:1100", "HGNC:1101"],
            target_type="gene_symbol",
        )

        assert result["HGNC:1100"] == "BRCA1"
        assert result["HGNC:1101"] == "BRCA2"

    def test_map_from_hgnc_to_ensembl(self, mapper):
        """Test mapping HGNC IDs to Ensembl IDs."""
        result = mapper.map_from_hgnc(
            ["HGNC:1100", "HGNC:1101"],
            target_type="ensembl_gene_id",
        )

        assert result["HGNC:1100"] == "ENSG00000012048"
        assert result["HGNC:1101"] == "ENSG00000139618"

    def test_map_ids_symbol_to_ensembl(self, mapper):
        """Test mapping between gene symbols and Ensembl IDs."""
        result = mapper.map_ids(
            ["BRCA1", "BRCA2"],
            source_type="gene_symbol",
            target_type="ensembl_gene_id",
        )

        assert result["BRCA1"] == "ENSG00000012048"
        assert result["BRCA2"] == "ENSG00000139618"

    def test_resolve_symbol_current(self, mapper):
        """Test resolving a current symbol."""
        result = mapper.resolve_symbol("BRCA1")
        assert result == "BRCA1"

    def test_resolve_symbol_alias(self, mapper):
        """Test resolving an alias symbol."""
        # BRCC1 is an alias for BRCA1
        result = mapper.resolve_symbol("BRCC1")
        assert result == "BRCA1"

        # FANCD1 is an alias for BRCA2
        result = mapper.resolve_symbol("FANCD1")
        assert result == "BRCA2"

    def test_resolve_symbol_unknown(self, mapper):
        """Test resolving an unknown symbol."""
        result = mapper.resolve_symbol("UNKNOWN_GENE")
        assert result is None

    def test_resolve_symbols_batch(self, mapper):
        """Test batch symbol resolution."""
        result = mapper.resolve_symbols(["BRCA1", "BRCC1", "FANCD1", "UNKNOWN"])

        assert result["BRCA1"] == "BRCA1"
        assert result["BRCC1"] == "BRCA1"
        assert result["FANCD1"] == "BRCA2"
        assert result["UNKNOWN"] is None

    def test_annotate_table_by_ensembl(self, mapper, hgnc_table):
        """Test annotating a table with HGNC data."""
        import hail as hl

        # Create a simple table with Ensembl IDs
        ht = hl.Table.parallelize([
            {"ensembl_id": "ENSG00000012048"},
            {"ensembl_id": "ENSG00000139618"},
            {"ensembl_id": "ENSG00000000000"},  # Unknown
        ]).key_by("ensembl_id")

        # Annotate
        result = mapper.annotate_table(
            ht,
            source_field="ensembl_id",
            source_type="ensembl_gene_id",
            fields_to_add=["hgnc_id", "gene_symbol"],
        )

        # Collect and check
        rows = result.collect()
        assert len(rows) == 3

        brca1_row = [r for r in rows if r.ensembl_id == "ENSG00000012048"][0]
        assert brca1_row.hgnc_hgnc_id == "HGNC:1100"
        assert brca1_row.hgnc_gene_symbol == "BRCA1"

    def test_get_gene_info(self, mapper):
        """Test getting full gene info."""
        info = mapper.get_gene_info("HGNC:1100")

        assert info is not None
        assert info["hgnc_id"] == "HGNC:1100"
        assert info["gene_symbol"] == "BRCA1"
        assert "BRCC1" in info["alias_symbols"]

    def test_get_gene_info_not_found(self, mapper):
        """Test getting gene info for non-existent gene."""
        info = mapper.get_gene_info("HGNC:9999999")
        assert info is None

    def test_get_genes_by_locus_group(self, mapper):
        """Test filtering genes by locus group."""
        protein_coding = mapper.get_genes_by_locus_group("protein-coding gene")

        # Test data has 3 protein-coding genes (BRCA1, BRCA2, GSTT1)
        assert len(protein_coding) == 3
        assert "HGNC:1100" in protein_coding
        assert "HGNC:1101" in protein_coding
        assert "HGNC:4641" in protein_coding

    def test_validate_ids(self, mapper):
        """Test ID validation."""
        result = mapper.validate_ids(
            ["BRCA1", "BRCA2", "UNKNOWN"],
            id_type="gene_symbol",
        )

        assert result["BRCA1"] is True
        assert result["BRCA2"] is True
        assert result["UNKNOWN"] is False

    def test_get_coverage_stats(self, mapper):
        """Test coverage statistics."""
        stats = mapper.get_coverage_stats()

        assert "total_genes" in stats
        assert stats["total_genes"] == 5  # 5 approved genes in test data

        assert "with_ensembl" in stats
        # BRCA1, BRCA2, GSTT1, FAM8A5P have Ensembl IDs
        assert stats["with_ensembl"] >= 4

        assert "protein_coding" in stats
        assert stats["protein_coding"] == 3
