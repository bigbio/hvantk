"""Utilities for loading and handling gene sets.

Provides a generic loader to replace the previous CHD-specific helper.
"""

from __future__ import annotations
from pathlib import Path
from typing import Iterable, Optional, Set, Union

__all__ = [
    "load_gene_set",
    "load_sample_chd_gene_set",
]


def load_gene_set(
    path: Optional[Union[str, Path]] = None,
    genes: Optional[Iterable[str]] = None,
    comment_prefix: str = "#",
    delimiter: Optional[str] = None,
    column: int = 0,
    strip_version: bool = False,
) -> Set[str]:
    """Load a gene set from a file and/or an iterable of gene symbols.

    Parameters
    ----------
    path : Optional[Union[str, Path]]
        Path to a text file containing gene identifiers. One gene per line (default) or
        delimited lines where the gene column is selected via `column`.
    genes : Optional[Iterable[str]]
        Additional genes to include (e.g. passed in code). Can be used without `path`.
    comment_prefix : str
        Lines beginning with this prefix are ignored.
    delimiter : Optional[str]
        Delimiter for splitting lines. If None, lines are treated as single tokens.
    column : int
        Column index (0-based) to extract when `delimiter` is provided.
    strip_version : bool
        If True, strip transcript / gene version suffix after '.' (e.g. ENSG0001.5 -> ENSG0001).

    Returns
    -------
    Set[str]
        Set of unique gene identifiers.

    Raises
    ------
    ValueError
        If neither path nor genes are provided.
    """
    gene_set: Set[str] = set()

    if path is None and genes is None:
        raise ValueError("Provide at least one of 'path' or 'genes'.")

    if path is not None:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"Gene set file not found: {p}")
        with p.open() as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith(comment_prefix):
                    continue
                if delimiter is not None:
                    parts = line.split(delimiter)
                    if column >= len(parts):
                        continue
                    token = parts[column].strip()
                else:
                    token = line
                if strip_version and "." in token:
                    token = token.split(".", 1)[0]
                if token:
                    gene_set.add(token)

    if genes is not None:
        for g in genes:
            token = g.strip()
            if strip_version and "." in token:
                token = token.split(".", 1)[0]
            if token:
                gene_set.add(token)

    return gene_set


def load_sample_chd_gene_set() -> Set[str]:
    """Return the legacy sample CHD-associated gene set (for backward compatibility)."""
    return {
        "GATA4",
        "NKX2-5",
        "TBX5",
        "NOTCH1",
        "CHD7",
        "TBX1",
        "MYH6",
        "ACTC1",
        "MYH7",
        "TNNT2",
    }
