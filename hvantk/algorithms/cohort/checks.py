"""Validate a cohort manifest against its on-disk table -- no Hail, no Click.

These checks used to be private to ``hvantk/tools/cohort/cohort_cli.py`` (design
decision D10 puts the cohort concept in ``hvantk/algorithms/cohort/``, but the CLI
module is the only place that could reach them). The consequence was concrete: the G1
acceptance gate (``local/planning/cohorts/g1_gate.py``) had to import three
underscore-private functions out of a CLI module, and any Python caller of
:func:`hvantk.algorithms.cohort.attach.attach` skipped this validation entirely.

Every function here raises plain ``ValueError`` -- this module must never import
``click``, or anything under ``hvantk.skills`` / ``hvantk.tools``. The CLI
(``hvantk/tools/cohort/cohort_cli.py``) keeps thin wrappers that convert
``ValueError`` into ``click.UsageError`` so the user-facing messages are unchanged.

Layering: stdlib only.
"""
from __future__ import annotations

import csv
import gzip
from pathlib import Path

_COMPRESSED_SUFFIXES = (".gz", ".bgz")


def _first_line(table_path: str) -> str:
    path = Path(table_path)
    if not path.exists():
        raise ValueError(f"cohort table not found: {table_path}")

    opener = gzip.open if path.suffix in _COMPRESSED_SUFFIXES else open
    with opener(path, "rt") as fh:
        line = fh.readline()
    if not line:
        raise ValueError(f"cohort table is empty: {table_path}")
    return line


def detect_delimiter(table_path: str) -> str:
    """Return the delimiter of a cohort table by sniffing its header line for a tab.

    Deliberately not derived from a fixed suffix whitelist (``.tsv``/``.txt`` vs
    everything else): that heuristic mis-delimits any tab-delimited file with an
    unlisted extension (``.tab``, no extension at all) or a compressed ``.tsv.gz``
    file, whose visible ``Path.suffix`` is ``.gz`` rather than ``.tsv``. Sniffing the
    header line for a literal tab is simpler and extension-agnostic. Gzip (``.gz`` /
    ``.bgz``) is decompressed transparently before sniffing.

    Used both by :func:`read_header` and by ``attach_cmd``'s ``hl.import_table`` call,
    so ``validate`` and ``attach`` can never disagree about the delimiter.
    """
    return "\t" if "\t" in _first_line(table_path) else ","


def is_compressed(table_path: str) -> bool:
    """Whether a cohort table is gzip-compressed, by the same suffix check
    :func:`_first_line` already applies.

    Public so a non-Hail reader of the table (``hvantk.algorithms.cohort.frame``'s
    pandas loader) can pick the same "is this gzip" answer as ``read_header`` and
    ``detect_delimiter`` without re-deriving the ``.gz``/``.bgz`` suffix list.
    """
    return Path(table_path).suffix in _COMPRESSED_SUFFIXES


def read_header(table_path: str) -> list[str]:
    """Return the column names of a delimited cohort table.

    Pure Python and header-only, so validation stays cheap and Hail-free. See
    :func:`detect_delimiter` for how the delimiter is chosen and how gzip is handled.
    """
    line = _first_line(table_path)
    delimiter = "\t" if "\t" in line else ","
    header = next(csv.reader([line], delimiter=delimiter))
    return [h.strip() for h in header]


def check_key_column_exists(manifest, header: list[str]) -> None:
    """The key column itself must exist in the table.

    ``key`` is the identifier SPACE (``gene_id`` / ``hgnc_id`` / ``symbol``) and
    ``key_column`` is the column that actually holds those values, and the two are
    commonly different (a real cohort's symbol column is often named ``gene``, not
    ``symbol``). Without this check, a wrong or missing ``key_column`` fails only deep
    inside ``hl.import_table`` with a raw Hail error that never mentions the manifest
    -- exactly the failure class every other declared column is already checked
    against.
    """
    if manifest.key_column not in header:
        raise ValueError(
            f"cohort table {manifest.table} has no column {manifest.key_column!r} "
            f"(key_column); key={manifest.key!r} names the identifier space, "
            "key_column names the column in the table that holds it"
        )


def check_declared_columns_exist(manifest, header: list[str]) -> None:
    """Every declared column must exist in the table.

    A missing declared column must never become a silently absent feature -- that is
    how a mapping bug turns into unexplained model degradation.
    """
    missing = [c for c in manifest.declared_columns() if c not in header]
    if missing:
        raise ValueError(
            f"cohort table {manifest.table} is missing declared column(s): "
            + ", ".join(missing)
        )
