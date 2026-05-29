"""GeneCatalogStreamer — ABC for gene metadata catalogs (HGNC, Ensembl, ...).

Concrete subclasses live in ``skills/<plugin>/streamers.py``. Algorithms
consume the ABC; callers wire concrete catalogs at construction time.
"""
from __future__ import annotations

from abc import ABC, abstractmethod

from hvantk.core.models import AnnotationTable


class GeneCatalogStreamer(ABC):
    """Per-DataModel operation contract for a gene metadata catalog."""

    def __init__(self, artifact: AnnotationTable) -> None:
        self._artifact = artifact

    @classmethod
    def from_path(cls, path: str) -> "GeneCatalogStreamer":
        """Construct a streamer from a persisted artifact path."""
        from hvantk.core.io import load
        return cls(load(path))

    @abstractmethod
    def is_canonical(self, symbol: str) -> bool:
        """Return True iff ``symbol`` is a currently approved symbol in this catalog."""

    @abstractmethod
    def resolve_alias(self, symbol: str) -> str | None:
        """Return canonical symbol if ``symbol`` is a known alias, else None."""

    @abstractmethod
    def expand_with_aliases(
        self, symbols: set[str]
    ) -> tuple[set[str], dict[str, str]]:
        """Expand ``symbols`` with known aliases.

        Returns
        -------
        (canonical_set, alias_to_canonical)
            canonical_set contains the input plus any resolved canonical
            forms; alias_to_canonical maps each input alias that was
            resolved to its canonical symbol.
        """

    @abstractmethod
    def map_ids(
        self,
        ids: list[str],
        source_type: str,
        target_type: str,
    ) -> dict[str, str | None]:
        """Map identifiers between ID types (e.g. gene_symbol -> hgnc_id).

        Used by gene-disease translation paths. Concrete catalogs (HGNC,
        Ensembl) implement the actual mapping.

        Parameters
        ----------
        ids : list of str
            Identifiers to translate.
        source_type : str
            ID type of inputs (e.g. "gene_symbol", "hgnc_id",
            "ensembl_id", "entrez_id", "uniprot_id").
        target_type : str
            ID type to translate to.

        Returns
        -------
        dict
            Mapping from each input id to its translated form (or None
            if untranslatable). Same keys as inputs.
        """

    def validate_symbols(
        self, symbols: set[str]
    ) -> tuple[set[str], dict[str, str], set[str]]:
        """Classify each symbol as canonical / alias-resolved / unrecognized.

        Returns
        -------
        (recognized, aliases_resolved, unrecognized)
        """
        recognized: set[str] = set()
        aliases_resolved: dict[str, str] = {}
        unrecognized: set[str] = set()
        for s in symbols:
            if self.is_canonical(s):
                recognized.add(s)
            else:
                canonical = self.resolve_alias(s)
                if canonical is not None:
                    aliases_resolved[s] = canonical
                    recognized.add(canonical)
                else:
                    unrecognized.add(s)
        return recognized, aliases_resolved, unrecognized
