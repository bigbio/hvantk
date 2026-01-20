"""
General-purpose Hail Table/MatrixTable utilities for field and schema handling.
"""

from typing import Union
import hail as hl


def get_row_fields(ht: Union[hl.Table, hl.MatrixTable]) -> set[str]:
    """
    Returns the set of row field names for a Hail Table or MatrixTable.
    Including key fields.

    Args:
        ht: A Hail Table or MatrixTable.

    Returns:
        A set of row field names as strings.
    """
    return set(ht.row.dtype.fields)


def get_col_fields(mt: hl.MatrixTable) -> set[str]:
    """
    Returns the set of column field names for a Hail MatrixTable.
    Including key fields.

    Args:
        mt: A Hail MatrixTable.

    Returns:
        A set of column field names as strings.
    """
    return set(mt.col.dtype.fields)


def get_entry_fields(mt: hl.MatrixTable) -> set[str]:
    """
    Returns the set of entry field names for a Hail MatrixTable.

    Args:
        mt: A Hail MatrixTable.

    Returns:
        A set of entry field names as strings.
    """
    return set(mt.entry.dtype.fields)
