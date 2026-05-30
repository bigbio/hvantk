"""VariantTableStreamer — ABC for variant table DataModels.

Per-variant rows keyed by (locus, alleles), arbitrary annotations.
Concrete subclasses live in ``skills/<plugin>/streamers.py``.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Iterable

import hail as hl

from hvantk.core.models import AnnotationTable


class VariantTableStreamer(ABC):
    """Per-DataModel operation contract for variant tables."""

    def __init__(self, artifact: AnnotationTable) -> None:
        self._artifact = artifact

    @classmethod
    def from_path(cls, path: str) -> "VariantTableStreamer":
        from hvantk.core.io import load
        return cls(load(path))

    @abstractmethod
    def to_hail(self) -> hl.Table:
        """Return the underlying Hail Table for downstream native operations."""

    @abstractmethod
    def filter_to_genes(self, genes: Iterable[str]) -> hl.Table:
        """Filter variants to those affecting any gene in ``genes``.

        'Affecting' is plugin-specific (ClinVar uses info.GENEINFO; dbNSFP a
        gene-symbol column).
        """

    @abstractmethod
    def filter_by_pathogenicity(self, labels: Iterable[str]) -> hl.Table:
        """Return variants whose pathogenicity label is in ``labels``.

        The label vocabulary is plugin-specific; the caller passes intended
        labels (e.g. ['Pathogenic', 'Likely_pathogenic']).
        """
