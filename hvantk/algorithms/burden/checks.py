"""Pure-Python contract checks for the gene-burden op (no Hail).

These validate the *names and vocabulary* of a burden request so a malformed call
fails fast with an actionable message before any Spark start-up, exactly like
``hvantk/algorithms/cohort/checks.py``. Data-dependent checks (bi-allelic rows,
binary arm) live in ``aggregate.assert_clean_mt`` because they need the MatrixTable.
"""
from __future__ import annotations

KEY_SPACES: tuple[str, ...] = ("gene_id", "hgnc_id", "symbol")
CARRIER_MODES: tuple[str, ...] = ("het", "hom", "chet", "homs_chet")


def check_key_space(key: str) -> None:
    if key not in KEY_SPACES:
        raise ValueError(f"unknown gene key space {key!r}; must be one of {KEY_SPACES}")


def check_carrier_mode(mode: str) -> None:
    if mode not in CARRIER_MODES:
        raise ValueError(
            f"unknown carrier mode {mode!r}; must be one of {CARRIER_MODES}"
        )


def check_required_fields(
    row_fields: list[str],
    entry_fields: list[str],
    col_fields: list[str],
    *,
    gene_col: str,
    route_col: str,
    arm_col: str,
) -> None:
    if gene_col not in row_fields:
        raise ValueError(
            f"gene column {gene_col!r} not found in MatrixTable row fields {row_fields}"
        )
    if route_col not in row_fields:
        raise ValueError(
            f"route column {route_col!r} not found in MatrixTable row fields {row_fields}"
        )
    if "GT" not in entry_fields:
        raise ValueError(
            f"genotype entry field 'GT' not found in MatrixTable entry fields {entry_fields}"
        )
    if arm_col not in col_fields:
        raise ValueError(
            f"case/control column {arm_col!r} not found in MatrixTable column fields {col_fields}"
        )
