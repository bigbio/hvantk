"""
Tests for table_utils module.
"""

import pytest
import hail as hl
from hvantk.core.utils.table_utils import (
    field_exists,
    get_row_fields,
    leaf_name,
    resolve_field,
    validate_fields,
)

# Mark as Hail-dependent
pytestmark = pytest.mark.hail


@pytest.fixture
def sample_table():
    """Create a sample Hail Table for testing."""
    # Create a simple table with various field types
    data = [
        {"id": 1, "name": "Alice", "age": 25, "score": 85.5},
        {"id": 2, "name": "Bob", "age": 30, "score": 92.0},
        {"id": 3, "name": "Charlie", "age": 35, "score": 78.3},
    ]

    ht = hl.Table.parallelize(data, key="id")
    return ht


@pytest.fixture
def sample_matrix_table():
    """Create a sample Hail MatrixTable for testing."""
    # Create a simple matrix table using a more straightforward approach
    n_rows, n_cols = 5, 3

    # Create the MatrixTable directly
    mt = hl.utils.range_matrix_table(n_rows, n_cols)

    # Add some row fields
    mt = mt.annotate_rows(
        gene_id=hl.str("gene_") + hl.str(mt.row_idx),
        chromosome=hl.str("chr") + hl.str((mt.row_idx % 3) + 1),
    )

    # Add some column fields
    mt = mt.annotate_cols(
        sample_id=hl.str("sample_") + hl.str(mt.col_idx),
        tissue=hl.str("tissue_") + hl.str(mt.col_idx % 2),
    )

    # Key by gene_id and sample_id
    mt = mt.key_rows_by("gene_id")
    mt = mt.key_cols_by("sample_id")

    return mt


@pytest.fixture
def empty_table():
    """Create an empty Hail Table for testing edge cases."""
    schema = hl.tstruct(id=hl.tint32, value=hl.tstr)
    return hl.Table.parallelize([], schema=schema, key="id")


def test_get_row_fields_table_basic(sample_table):
    """Test get_row_fields with a basic Hail Table."""
    fields = get_row_fields(sample_table)

    # Check that it returns a set
    assert isinstance(fields, set)

    # Check that all expected fields are present
    expected_fields = {"id", "name", "age", "score"}
    assert fields == expected_fields


def test_get_row_fields_matrix_table_basic(sample_matrix_table):
    """Test get_row_fields with a basic Hail MatrixTable."""
    fields = get_row_fields(sample_matrix_table)

    # Check that it returns a set
    assert isinstance(fields, set)

    # Check that all expected row fields are present (including row_idx from range_matrix_table)
    expected_fields = {"gene_id", "chromosome", "row_idx"}
    assert fields == expected_fields


def test_get_row_fields_empty_table(empty_table):
    """Test get_row_fields with an empty table."""
    fields = get_row_fields(empty_table)

    # Check that it returns a set
    assert isinstance(fields, set)

    # Check that expected fields are present even for empty table
    expected_fields = {"id", "value"}
    assert fields == expected_fields


def test_get_row_fields_table_with_key_fields():
    """Test that key fields are included in the result."""
    # Create a table with multiple key fields
    data = [
        {"chr": "chr1", "pos": 1000, "ref": "A", "alt": "T", "qual": 30.0},
        {"chr": "chr1", "pos": 2000, "ref": "G", "alt": "C", "qual": 25.5},
        {"chr": "chr2", "pos": 1500, "ref": "C", "alt": "G", "qual": 40.2},
    ]

    ht = hl.Table.parallelize(data, key=["chr", "pos"])
    fields = get_row_fields(ht)

    # Check that key fields are included
    expected_fields = {"chr", "pos", "ref", "alt", "qual"}
    assert fields == expected_fields


def test_get_row_fields_table_with_complex_types():
    """Test get_row_fields with complex field types."""
    # Create a table with complex types
    data = [
        {
            "id": 1,
            "info": {"AC": 2, "AF": 0.5},
            "samples": ["sample1", "sample2"],
            "coordinates": hl.Struct(x=1.0, y=2.0),
        },
        {
            "id": 2,
            "info": {"AC": 1, "AF": 0.25},
            "samples": ["sample3"],
            "coordinates": hl.Struct(x=3.0, y=4.0),
        },
    ]

    ht = hl.Table.parallelize(data, key="id")
    fields = get_row_fields(ht)

    # Check that all fields are present regardless of type complexity
    expected_fields = {"id", "info", "samples", "coordinates"}
    assert fields == expected_fields


def test_get_row_fields_matrix_table_with_complex_row_fields():
    """Test get_row_fields with MatrixTable having complex row fields."""
    # Create a simple MatrixTable and add complex row annotations
    mt = hl.utils.range_matrix_table(3, 2)

    # Add complex row fields
    mt = mt.annotate_rows(
        variant_id=hl.str("variant_") + hl.str(mt.row_idx),
        locus=hl.locus("chr1", 1000 + mt.row_idx * 100),
        alleles=hl.array(["A", "T"]),
        info=hl.struct(AC=5, AN=10),
        filters=hl.set(["PASS"]),
    )

    # Key by variant_id
    mt = mt.key_rows_by("variant_id")

    fields = get_row_fields(mt)

    # Check that all row fields are present (including row_idx from range_matrix_table)
    expected_fields = {"variant_id", "locus", "alleles", "info", "filters", "row_idx"}
    assert fields == expected_fields


def test_get_row_fields_return_type():
    """Test that get_row_fields always returns a set of strings."""
    # Create a simple table
    data = [{"x": 1, "y": "test"}]
    ht = hl.Table.parallelize(data, key="x")

    fields = get_row_fields(ht)

    # Check return type
    assert isinstance(fields, set)

    # Check that all elements are strings
    assert all(isinstance(field, str) for field in fields)


def test_get_row_fields_consistency(sample_table):
    """Test that get_row_fields returns consistent results for the same table."""
    # Call multiple times on the same table
    fields1 = get_row_fields(sample_table)
    fields2 = get_row_fields(sample_table)
    fields3 = get_row_fields(sample_table)

    # Results should be identical
    assert fields1 == fields2 == fields3


def test_get_row_fields_different_table_types():
    """Test get_row_fields with different table configurations."""
    # Table with no key
    data1 = [{"field1": "a", "field2": 1}]
    ht1 = hl.Table.parallelize(data1)  # No key specified
    fields1 = get_row_fields(ht1)
    expected1 = {"field1", "field2"}
    assert fields1 == expected1

    # Table with single key
    data2 = [{"id": 1, "value": "test"}]
    ht2 = hl.Table.parallelize(data2, key="id")
    fields2 = get_row_fields(ht2)
    expected2 = {"id", "value"}
    assert fields2 == expected2

    # Table with compound key
    data3 = [{"chr": "1", "pos": 100, "data": "example"}]
    ht3 = hl.Table.parallelize(data3, key=["chr", "pos"])
    fields3 = get_row_fields(ht3)
    expected3 = {"chr", "pos", "data"}
    assert fields3 == expected3


# Tests for get_col_fields function
def test_get_col_fields_basic(sample_matrix_table):
    """Test get_col_fields with a basic MatrixTable."""
    from hvantk.core.utils.table_utils import get_col_fields

    fields = get_col_fields(sample_matrix_table)

    # Check that it returns a set
    assert isinstance(fields, set)

    # Check that all expected column fields are present
    expected_fields = {"sample_id", "tissue", "col_idx"}
    assert fields == expected_fields


def test_get_col_fields_with_complex_types():
    """Test get_col_fields with complex column field types."""
    from hvantk.core.utils.table_utils import get_col_fields

    # Create a MatrixTable with complex column fields
    mt = hl.utils.range_matrix_table(2, 3)

    mt = mt.annotate_cols(
        sample_info=hl.struct(
            id=hl.str("sample_") + hl.str(mt.col_idx),
            metadata=hl.dict({"batch": "batch1", "tissue": "brain"}),
        ),
        phenotypes=hl.array([1.5, 2.3, 0.8]),
        flags=hl.set(["QC_PASS", "HIGH_QUALITY"]),
    )

    fields = get_col_fields(mt)

    # Check that all fields are present
    expected_fields = {"col_idx", "sample_info", "phenotypes", "flags"}
    assert fields == expected_fields


def test_get_entry_fields_basic():
    """Test get_entry_fields with a basic MatrixTable."""
    from hvantk.core.utils.table_utils import get_entry_fields

    # Create a MatrixTable with entry fields
    mt = hl.utils.range_matrix_table(2, 2)

    # Add some entry fields
    mt = mt.annotate_entries(
        genotype=hl.call(0, 1),
        depth=mt.row_idx + mt.col_idx + 10,
        quality=hl.float64(30.5),
    )

    fields = get_entry_fields(mt)

    # Check that it returns a set
    assert isinstance(fields, set)

    # Check that all expected entry fields are present
    expected_fields = {"genotype", "depth", "quality"}
    assert fields == expected_fields


def test_get_entry_fields_empty():
    """Test get_entry_fields with a MatrixTable with no entry fields."""
    from hvantk.core.utils.table_utils import get_entry_fields

    # Create a basic MatrixTable with no additional entry fields
    mt = hl.utils.range_matrix_table(2, 2)

    fields = get_entry_fields(mt)

    # Check that it returns an empty set
    assert isinstance(fields, set)
    assert len(fields) == 0


# ---------------------------------------------------------------------------
# Tests for resolve_field / field_exists / validate_fields
# ---------------------------------------------------------------------------


@pytest.fixture
def table_with_struct():
    """Table with a nested struct field."""
    data = [
        {"id": 1, "phe": hl.Struct(is_case=True, age=25)},
        {"id": 2, "phe": hl.Struct(is_case=False, age=30)},
    ]
    return hl.Table.parallelize(data, key="id")


@pytest.fixture
def table_with_dotted_field():
    """Table with a flat field whose name contains a dot (as from flatten())."""
    ht = hl.Table.parallelize(
        [{"id": 1, "phe": hl.Struct(is_case=True, age=25)}],
        key="id",
    )
    # flatten() converts struct fields into dot-delimited flat fields
    return ht.flatten()


def test_resolve_field_flat(sample_table):
    """resolve_field finds a plain top-level field."""
    expr = resolve_field(sample_table, "name")
    assert expr.dtype == hl.tstr


def test_resolve_field_flat_dotted(table_with_dotted_field):
    """resolve_field finds a flat field whose name contains a dot."""
    # After flatten(), the field is literally named "phe.is_case"
    expr = resolve_field(table_with_dotted_field, "phe.is_case")
    assert expr.dtype == hl.tbool


def test_resolve_field_struct_navigation(table_with_struct):
    """resolve_field navigates into a struct when no flat match exists."""
    expr = resolve_field(table_with_struct, "phe.is_case")
    assert expr.dtype == hl.tbool


def test_resolve_field_struct_navigation_deep():
    """resolve_field navigates multiple struct levels."""
    data = [
        {"id": 1, "a": hl.Struct(b=hl.Struct(c=42))},
    ]
    ht = hl.Table.parallelize(data, key="id")
    expr = resolve_field(ht, "a.b.c")
    assert hl.eval(expr.collect()[0]) == 42


def test_resolve_field_missing_raises(sample_table):
    """resolve_field raises LookupError with available fields listed."""
    with pytest.raises(LookupError, match="not found"):
        resolve_field(sample_table, "nonexistent")


def test_resolve_field_missing_dotted_raises(sample_table):
    """resolve_field raises when both literal and struct navigation fail."""
    with pytest.raises(LookupError, match="not found"):
        resolve_field(sample_table, "no.such.field")


def test_resolve_field_on_struct_expr(table_with_struct):
    """resolve_field works directly on a StructExpression (e.g., mt.row)."""
    expr = resolve_field(table_with_struct.row, "phe.is_case")
    assert expr.dtype == hl.tbool


def test_field_exists_true(sample_table):
    """field_exists returns True for an existing field."""
    assert field_exists(sample_table, "name") is True


def test_field_exists_struct(table_with_struct):
    """field_exists returns True for a nested struct path."""
    assert field_exists(table_with_struct, "phe.is_case") is True


def test_field_exists_false(sample_table):
    """field_exists returns False for a missing field."""
    assert field_exists(sample_table, "nonexistent") is False


def test_validate_fields_all_valid(sample_table):
    """validate_fields returns empty list when all fields exist."""
    errors = validate_fields(sample_table, ["name", "age", "score"])
    assert errors == []


def test_validate_fields_some_missing(sample_table):
    """validate_fields returns errors for missing fields only."""
    errors = validate_fields(
        sample_table,
        ["name", "missing1", "age", "missing2"],
        context="test table",
    )
    assert len(errors) == 2
    assert "missing1" in errors[0]
    assert "test table" in errors[0]
    assert "missing2" in errors[1]


def test_validate_fields_struct_path(table_with_struct):
    """validate_fields accepts dotted struct paths."""
    errors = validate_fields(table_with_struct, ["phe.is_case", "phe.age"])
    assert errors == []


# ---------------------------------------------------------------------------
# Tests for leaf_name (no Hail needed)
# ---------------------------------------------------------------------------


class TestLeafName:
    """Tests for leaf_name — pure string utility, no Hail required."""

    def test_dotted_path(self):
        assert leaf_name("phe.is_case") == "is_case"

    def test_multi_level(self):
        assert leaf_name("a.b.c") == "c"

    def test_no_dots(self):
        assert leaf_name("status") == "status"

    def test_empty_string(self):
        assert leaf_name("") == ""
