"""Tests for global metadata annotation on built Hail Tables and MatrixTables."""

import pytest

pytestmark = [pytest.mark.hail]


def test_build_table_metadata_fields(tmp_path):
    """Test that build_table_metadata produces correct struct fields."""
    import hail as hl
    from hvantk.core.models.metadata import build_table_metadata

    input_path = str(tmp_path / "input.tsv")
    with open(input_path, "w") as f:
        f.write("gene\tscore\n")
        f.write("BRCA1\t0.99\n")

    ht = hl.import_table(input_path, types={"score": hl.tfloat64}).key_by("gene")
    metadata = build_table_metadata("TestSource", input_path, ht)

    # Annotate and evaluate globals
    ht = ht.annotate_globals(hvantk_metadata=metadata)
    meta = hl.eval(ht.hvantk_metadata)

    assert meta.source_name == "TestSource"
    assert meta.raw_input_path == input_path
    assert "gene" in meta.row_schema
    assert "score" in meta.row_schema
    assert meta.key_fields == ["gene"]
    assert "gene" in meta.key_schema
    assert meta.n_fields >= 2
    assert meta.build_date is not None
    assert meta.hvantk_version is not None
    assert meta.reference_genome == "NA"  # no locus field


def test_build_table_metadata_with_locus():
    """Test metadata captures reference genome when locus is present."""
    import hail as hl
    from hvantk.core.models.metadata import build_table_metadata

    ht = hl.utils.range_table(5)
    ht = ht.annotate(locus=hl.locus("chr1", ht.idx + 1, reference_genome="GRCh38"))
    ht = ht.key_by("locus")

    metadata = build_table_metadata("LocusTest", "/fake/path.vcf", ht)
    ht = ht.annotate_globals(hvantk_metadata=metadata)
    meta = hl.eval(ht.hvantk_metadata)

    assert meta.reference_genome == "GRCh38"


def test_build_matrix_metadata_fields(tmp_path):
    """Test that build_matrix_metadata produces correct struct fields."""
    import hail as hl
    from hvantk.core.models.metadata import build_matrix_metadata

    mt = hl.utils.range_matrix_table(2, 3)
    mt = mt.annotate_rows(rid=mt.row_idx).key_rows_by("rid")
    mt = mt.annotate_cols(sample=hl.str(mt.col_idx)).key_cols_by("sample")

    input_path = str(tmp_path / "matrix.tsv")
    metadata = build_matrix_metadata("ExpressionAtlas", input_path, mt)
    mt = mt.annotate_globals(hvantk_metadata=metadata)
    meta = hl.eval(mt.hvantk_metadata)

    assert meta.source_name == "ExpressionAtlas"
    assert meta.source_description != ""
    assert meta.raw_input_path == input_path
    assert "rid" in meta.row_schema
    assert "sample" in meta.col_schema
    assert meta.key_fields == ["rid"]
    assert meta.n_row_fields >= 1
    assert meta.n_col_fields >= 1
    assert meta.n_cols is None


def test_get_hvantk_version():
    """Test version retrieval function."""
    from hvantk.core.models.metadata import _get_hvantk_version

    version = _get_hvantk_version()
    assert isinstance(version, str)
    assert len(version) > 0
