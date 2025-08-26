"""
General-purpose Hail Table/MatrixTable utilities for field and schema handling.
"""

from typing import Union
import hail as hl

def get_row_fields(ht: Union[hl.Table, hl.MatrixTable]) -> set[str]:
    """
    Returns the set of row field names for a Hail Table or MatrixTable.

    Args:
        ht: A Hail Table or MatrixTable.

    Returns:
        A set of row field names as strings.
    """
    return set(ht.row.dtype.fields.keys())
