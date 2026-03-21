"""Build metadata structs for annotating Hail Tables and MatrixTables globals."""

from datetime import datetime

import hail as hl


def _get_hvantk_version() -> str:
    """Get hvantk version from package metadata."""
    try:
        from importlib.metadata import version

        return version("hvantk")
    except Exception:
        return "unknown"


def build_table_metadata(
    source_name: str,
    input_path: str,
    ht: hl.Table,
) -> hl.struct:
    """Build an hvantk_metadata struct for a Hail Table.

    Parameters
    ----------
    source_name : str
        Human-readable name of the data source.
    input_path : str
        Path to the raw input file.
    ht : hl.Table
        The Hail Table (used to extract schema info).

    Returns
    -------
    hl.struct
        Metadata struct suitable for ``ht.annotate_globals()``.
    """
    return hl.struct(
        hvantk_version=_get_hvantk_version(),
        source_name=source_name,
        raw_input_path=input_path,
        build_date=datetime.now().isoformat(),
        reference_genome=str(ht.locus.dtype.reference_genome)
        if "locus" in ht.row
        else "NA",
        row_schema=str(ht.row.dtype),
        key_schema=str(ht.key.dtype),
        key_fields=list(ht.key),
        n_fields=len(ht.row),
    )


def build_matrix_metadata(
    source_name: str,
    input_path: str,
    mt: hl.MatrixTable,
) -> hl.struct:
    """Build an hvantk_metadata struct for a Hail MatrixTable.

    Parameters
    ----------
    source_name : str
        Human-readable name of the data source.
    input_path : str
        Path to the raw input file.
    mt : hl.MatrixTable
        The Hail MatrixTable (used to extract schema info).

    Returns
    -------
    hl.struct
        Metadata struct suitable for ``mt.annotate_globals()``.
    """
    return hl.struct(
        hvantk_version=_get_hvantk_version(),
        source_name=source_name,
        raw_input_path=input_path,
        build_date=datetime.now().isoformat(),
        row_schema=str(mt.row.dtype),
        col_schema=str(mt.col.dtype),
        entry_schema=str(mt.entry.dtype),
        key_schema=str(mt.row_key.dtype),
        key_fields=list(mt.row_key),
        n_row_fields=len(mt.row),
        n_col_fields=len(mt.col),
    )
