"""Genome / contig related utility functions.

Currently includes helpers for constructing contig recoding dictionaries
for converting plain numeric (and sex) chromosome labels (e.g. '1', '2', 'X')
into reference-style contigs (e.g. 'chr1', 'chr2', 'chrX') when importing
VCF / BED data into Hail.
"""
from __future__ import annotations
from typing import Iterable, Tuple, Dict

__all__ = [
    "contig_recoding",
]

def contig_recoding(
    prefix: str = "chr",
    autosomes: Iterable[int] = range(1, 23),
    sex_chromosomes: Iterable[str] = ("X", "Y"),
    include_mt: bool = False,
    source_mt_names: Tuple[str, ...] = ("MT", "M"),
    target_mt_name: str | None = None,
) -> Dict[str, str]:
    """Return a contig recoding dict for Hail import_* functions.

    Parameters
    ----------
    prefix : str
        Prefix to prepend to chromosome names, default 'chr'. Use '' for no prefix.
    autosomes : Iterable[int]
        Collection of autosome numbers to include (default 1..22).
    sex_chromosomes : Iterable[str]
        Sex chromosome labels to include.
    include_mt : bool
        Whether to include mitochondrial chromosome mapping.
    source_mt_names : Tuple[str, ...]
        Possible source contig labels for mitochondrial chromosome in the input data.
    target_mt_name : str | None
        Target (post-recode) mitochondrial contig suffix (without prefix). If None and prefix=='chr', uses 'M';
        else uses first element of `source_mt_names`.

    Returns
    -------
    Dict[str, str]
        Mapping suitable for the `contig_recoding` parameter of Hail import functions.

    Examples
    --------
    >>> contig_recoding()
    {'1': 'chr1', '2': 'chr2', ..., '22': 'chr22', 'X': 'chrX', 'Y': 'chrY'}
    >>> contig_recoding(include_mt=True)['MT']
    'chrM'
    """
    recode: Dict[str, str] = {}
    # Autosome mapping
    for a in autosomes:
        recode[str(a)] = f"{prefix}{a}" if prefix else str(a)
    # Sex chromosomes
    for sc in sex_chromosomes:
        recode[str(sc)] = f"{prefix}{sc}" if prefix else str(sc)
    # Mitochondrial chromosome(s)
    if include_mt:
        if target_mt_name is None:
            if prefix == "chr":
                target_mt_name = "M"
            else:
                target_mt_name = source_mt_names[0]
        target = f"{prefix}{target_mt_name}" if prefix else target_mt_name
        for src in source_mt_names:
            recode[src] = target
    return recode
