"""
General-purpose Hail Table/MatrixTable utilities for field and schema handling.
"""

#TODO: All table utils should go here

def get_row_fields(ht) -> set[str]:
    """
    Returns the set of row field names for a Hail Table or MatrixTable.

    Args:
        ht: A Hail Table or MatrixTable.

    Returns:
        A set of row field names as strings.
    """
    return set(ht.row.dtype.fields)
