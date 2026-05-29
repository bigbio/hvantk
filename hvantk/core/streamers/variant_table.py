"""VariantTableStreamer — ABC for variant table DataModels.

Per-variant rows keyed by (locus, alleles), arbitrary annotations.

Concrete subclasses live in ``skills/<plugin>/streamer.py``. This stub
is filled with the full abstract surface in a later PR, when the first
concrete consumer (ClinVarVariantTableStreamer) is built.
"""
from __future__ import annotations

from abc import ABC, abstractmethod

import hail as hl

from hvantk.core.models import AnnotationTable


class VariantTableStreamer(ABC):
    """Per-DataModel operation contract for variant tables (stub)."""

    def __init__(self, artifact: AnnotationTable) -> None:
        self._artifact = artifact

    @classmethod
    def from_path(cls, path: str) -> "VariantTableStreamer":
        """Construct a streamer from a persisted artifact path."""
        from hvantk.core.io import load
        return cls(load(path))

    @abstractmethod
    def to_hail(self) -> hl.Table:
        """Return the underlying Hail Table.

        Stub method retained so subclasses must provide a Hail-level handle.
        Additional abstract methods land in a later PR.
        """
