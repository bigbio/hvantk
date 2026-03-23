"""
General-purpose Hail Table/MatrixTable utilities for field and schema handling.
"""

from __future__ import annotations

import logging
import re
from typing import Any, Dict, List, Optional, Sequence, Union

import hail as hl

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Field introspection
# ---------------------------------------------------------------------------


def _normalize_field_name(name: str) -> str:
    """Normalize a field name for case/separator-insensitive comparison.

    Strips all non-alphanumeric characters and lowercases, so that
    ``"Gene Symbol"``, ``"GENE_SYMBOL"``, and ``"gene_symbol"`` all
    produce the same key ``"genesymbol"``.
    """
    return re.sub(r"[^a-z0-9]", "", name.lower())


def build_rename_map(
    field_map: Dict[str, str],
    actual_fields: set[str],
) -> Dict[str, str]:
    """Build a rename map from *field_map* to match *actual_fields*.

    Performs **case-insensitive, separator-insensitive** matching so that
    rename maps written for one header convention (e.g. title-case with
    spaces) also work for other conventions (e.g. UPPER_SNAKE_CASE).

    Parameters
    ----------
    field_map : dict
        Canonical rename map ``{raw_header: target_name}``.
    actual_fields : set of str
        Field names actually present in the table.

    Returns
    -------
    dict
        ``{actual_field: target_name}`` for every matched field.
        Fields already matching a target name are skipped.
    """
    # Index actual fields by their normalized form
    norm_to_actual: Dict[str, str] = {}
    for f in actual_fields:
        norm_to_actual[_normalize_field_name(f)] = f

    target_names = set(field_map.values())
    rename = {}
    for canonical_raw, target in field_map.items():
        norm = _normalize_field_name(canonical_raw)
        actual = norm_to_actual.get(norm)
        if actual is not None and actual != target:
            # Skip if actual field name already equals the target or
            # if the target name is already taken by another field
            if actual not in target_names:
                rename[actual] = target
        # Also try matching by the target name's normalized form
        # (handles cases where actual headers already use the target naming)
        if actual is None:
            norm_target = _normalize_field_name(target)
            actual_by_target = norm_to_actual.get(norm_target)
            if actual_by_target is not None and actual_by_target != target:
                if actual_by_target not in target_names:
                    rename[actual_by_target] = target

    if rename:
        logger.debug("Header rename map: %s", rename)
    return rename


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


# ---------------------------------------------------------------------------
# String → Boolean coercion
# ---------------------------------------------------------------------------

_TRUTHY_VALUES = hl.set({"yes", "y", "true", "1"})


def str_to_bool(expr: hl.StringExpression) -> hl.BooleanExpression:
    """Coerce a string field to boolean using common truthy values.

    Recognises (case-insensitive): ``"yes"``, ``"y"``, ``"true"``, ``"1"``.
    Missing or any other value maps to ``False``.

    Parameters
    ----------
    expr : hl.StringExpression
        The string field to convert.

    Returns
    -------
    hl.BooleanExpression
    """
    return hl.if_else(
        hl.is_defined(expr) & _TRUTHY_VALUES.contains(expr.lower()),
        True,
        False,
    )


# ---------------------------------------------------------------------------
# Field resolution — handles both flat dotted names and nested structs
# ---------------------------------------------------------------------------


def resolve_field(obj: Any, field_path: str) -> Any:
    """Resolve a field path on a Hail Table, MatrixTable, or StructExpression.

    Supports two naming styles transparently:

    1. **Flat dotted names** — a field literally named ``"phe.is_case"``
       (common after ``Table.flatten()`` or in user-annotated schemas).
       Accessed via ``obj["phe.is_case"]``.

    2. **Nested struct paths** — a struct ``phe`` with sub-field ``is_case``.
       Accessed via ``obj["phe"]["is_case"]``.

    The function tries a literal lookup first.  If the field is not found and
    ``field_path`` contains dots, it falls back to chained struct navigation.

    Parameters
    ----------
    obj
        A Hail :class:`~hail.Table`, :class:`~hail.MatrixTable`, or
        :class:`~hail.expr.StructExpression`.
    field_path : str
        Field name — may contain dots for either flat or nested access.

    Returns
    -------
    hail.expr.Expression
        The resolved Hail expression.

    Raises
    ------
    LookupError
        If the field cannot be found via either strategy.
    """
    # --- Strategy 1: literal lookup (handles flat dotted names) ---
    try:
        return obj[field_path]
    except (LookupError, KeyError):
        pass

    # --- Strategy 2: struct navigation (handles nested structs) ---
    if "." in field_path:
        parts = field_path.split(".")
        try:
            expr = obj
            for part in parts:
                expr = expr[part]
            logger.debug(
                "Resolved '%s' via struct navigation (not a flat field)",
                field_path,
            )
            return expr
        except (LookupError, KeyError):
            pass

    # --- Both strategies failed — build a helpful error message ---
    available = _available_fields(obj)
    available_str = ", ".join(sorted(available)) if available else "(none)"
    raise LookupError(
        f"Field '{field_path}' not found. "
        f"Available top-level fields: {available_str}"
    )


def field_exists(obj: Any, field_path: str) -> bool:
    """Check whether *field_path* can be resolved on *obj*.

    Uses :func:`resolve_field` internally, so it handles both flat dotted
    names and nested struct paths.

    Parameters
    ----------
    obj
        A Hail Table, MatrixTable, or StructExpression.
    field_path : str
        Field name — may contain dots.

    Returns
    -------
    bool
    """
    try:
        resolve_field(obj, field_path)
        return True
    except (LookupError, KeyError):
        return False


def leaf_name(field_path: str) -> str:
    """Return the last component of a dot-delimited field path.

    Examples
    --------
    >>> leaf_name("phe.is_case")
    'is_case'
    >>> leaf_name("status")
    'status'
    """
    return field_path.rsplit(".", 1)[-1]


def validate_fields(
    obj: Any,
    field_paths: Sequence[str],
    *,
    context: str = "",
) -> List[str]:
    """Validate that all *field_paths* are resolvable on *obj*.

    Parameters
    ----------
    obj
        A Hail Table, MatrixTable, or StructExpression.
    field_paths
        Field names to validate.
    context : str, optional
        Label for error messages (e.g. ``"cohort MT row fields"``).

    Returns
    -------
    List[str]
        Error messages for fields that could not be resolved.
        Empty list means all fields are valid.
    """
    errors: List[str] = []
    available: Optional[set] = None
    for fp in field_paths:
        if not field_exists(obj, fp):
            if available is None:
                available = _available_fields(obj)
            ctx = f" in {context}" if context else ""
            errors.append(
                f"Field '{fp}' not found{ctx}. "
                f"Available: {', '.join(sorted(available)) if available else '(none)'}"
            )
    return errors


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _available_fields(obj: Any) -> set:
    """Return the set of top-level field names for *obj*."""
    # StructExpression / row / col / entry
    if hasattr(obj, "dtype") and hasattr(obj.dtype, "fields"):
        return set(obj.dtype.fields)
    # Table
    if isinstance(obj, hl.Table):
        return set(obj.row.dtype.fields)
    # MatrixTable — return row + col + entry for maximum helpfulness
    if isinstance(obj, hl.MatrixTable):
        return (
            set(obj.row.dtype.fields)
            | set(obj.col.dtype.fields)
            | set(obj.entry.dtype.fields)
        )
    return set()
