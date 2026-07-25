"""Identifier reconciliation onto the gene spine.

Every annotation source arrives keyed on something other than the spine's ``gene_id``:
HGNC IDs (ClinGen, GenCC), gene symbols (expression atlases), or genomic coordinates.
Reconciliation happens here and nowhere else, and it always produces a
:class:`MappingReport` so loss is measured rather than silent.

The failure this exists to prevent: the previous per-gene matrix merged on raw symbols
with no alias resolution, so every previous-symbol or alias mismatch dropped a gene and
nothing reported it.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Iterable, Optional

logger = logging.getLogger(__name__)

MAX_NAMED_UNMAPPED = 20


class MappingRateError(Exception):
    """A source mapped below its declared minimum rate."""


@dataclass(frozen=True)
class MappingReport:
    """How much of a source reached the spine, and what did not."""

    source: str
    key_type: str
    n_in: int
    n_mapped: int
    n_unmapped: int
    n_ambiguous: int
    unmapped: tuple[str, ...]

    @property
    def rate(self) -> float:
        """Mapped fraction. An empty input maps completely, vacuously."""
        return self.n_mapped / self.n_in if self.n_in else 1.0

    def summary(self) -> str:
        return (
            f"{self.source} [{self.key_type}]: {self.n_mapped}/{self.n_in} mapped "
            f"({self.rate:.2%}), {self.n_unmapped} unmapped, "
            f"{self.n_ambiguous} ambiguous"
        )


def enforce_rate(report: MappingReport, min_rate: float) -> None:
    """Raise if ``report`` fell below ``min_rate``.

    A source that maps poorly must fail loudly. Contributing a column of NaNs instead is
    how mapping bugs become invisible model degradation.
    """
    if report.rate < min_rate:
        examples = ", ".join(report.unmapped[:10])
        raise MappingRateError(
            f"{report.source} mapped {report.rate:.2f} of its identifiers, below the "
            f"required {min_rate:.2f}. Unmapped examples: {examples}"
        )


class GeneIdMapper:
    """Map source identifiers onto spine ``gene_id`` values.

    Parameters
    ----------
    hgnc : HGNCGeneCatalogStreamer
        Built from the ``hgnc:lookup`` artifact.
    spine_gene_ids : set of str
        The spine's ``gene_id`` values. An identifier that resolves to an Ensembl ID
        outside the spine counts as unmapped, not mapped -- otherwise a join would
        silently drop it later, which is the loss this class exists to surface.
    """

    def __init__(self, hgnc, spine_gene_ids: Iterable[str]) -> None:
        self._hgnc = hgnc
        self._spine = set(spine_gene_ids)

    def from_hgnc_ids(
        self, hgnc_ids: Iterable[str], *, source: str
    ) -> tuple[dict[str, Optional[str]], MappingReport]:
        ids = list(hgnc_ids)
        to_ensembl = self._hgnc.map_from_hgnc(ids, "ensembl_gene_id")
        mapping = {i: self._keep_if_on_spine(to_ensembl.get(i)) for i in ids}
        return mapping, self._report(source, "hgnc_id", mapping)

    def from_gene_ids(
        self, gene_ids: Iterable[str], *, source: str
    ) -> tuple[dict[str, Optional[str]], MappingReport]:
        """Resolve source rows already keyed on Ensembl ``gene_id``.

        No HGNC lookup: the identifiers are already in the spine's key space. A
        ``gene_id`` the spine does not contain (a different Ensembl release, or a
        non-protein-coding gene) counts as unmapped, exactly as for the other key types.
        """
        ids = list(gene_ids)
        mapping = {g: self._keep_if_on_spine(g) for g in ids}
        return mapping, self._report(source, "gene_id", mapping)

    def from_symbols(
        self, symbols: Iterable[str], *, source: str
    ) -> tuple[dict[str, Optional[str]], MappingReport]:
        syms = list(symbols)
        canonical = {s: self._hgnc.resolve_to_canonical(s) for s in syms}
        to_hgnc = self._hgnc.map_to_hgnc(
            sorted({c for c in canonical.values() if c}), "gene_symbol"
        )
        hgnc_ids = sorted({h for h in to_hgnc.values() if h})
        to_ensembl = self._hgnc.map_from_hgnc(hgnc_ids, "ensembl_gene_id")

        mapping: dict[str, Optional[str]] = {}
        for sym in syms:
            hgnc_id = to_hgnc.get(canonical.get(sym))
            gene_id = to_ensembl.get(hgnc_id) if hgnc_id else None
            mapping[sym] = self._keep_if_on_spine(gene_id)
        return mapping, self._report(source, "symbol", mapping)

    def from_uniprot_ids(
        self, uniprot_ids: Iterable[str], *, source: str
    ) -> tuple[dict[str, Optional[str]], MappingReport]:
        """Resolve rows keyed on UniProt accession (proteomic / interactome sources).

        Same two-hop chain as ``from_symbols`` (accession -> hgnc_id -> ensembl_gene_id)
        but with no alias step: UniProt accessions are stable, so there is nothing to
        canonicalise. An accession HGNC does not carry -- an isoform-specific or
        non-human entry -- counts as unmapped.
        """
        ids = list(uniprot_ids)
        to_hgnc = self._hgnc.map_to_hgnc(sorted(set(ids)), "uniprot_id")
        hgnc_ids = sorted({h for h in to_hgnc.values() if h})
        to_ensembl = self._hgnc.map_from_hgnc(hgnc_ids, "ensembl_gene_id")

        mapping: dict[str, Optional[str]] = {}
        for acc in ids:
            hgnc_id = to_hgnc.get(acc)
            gene_id = to_ensembl.get(hgnc_id) if hgnc_id else None
            mapping[acc] = self._keep_if_on_spine(gene_id)
        return mapping, self._report(source, "uniprot_id", mapping)

    def _keep_if_on_spine(self, gene_id: Optional[str]) -> Optional[str]:
        return gene_id if gene_id in self._spine else None

    def _report(
        self, source: str, key_type: str, mapping: dict[str, Optional[str]]
    ) -> MappingReport:
        unmapped = tuple(sorted(k for k, v in mapping.items() if v is None))
        targets = [v for v in mapping.values() if v is not None]
        n_ambiguous = len(targets) - len(set(targets))
        report = MappingReport(
            source=source,
            key_type=key_type,
            n_in=len(mapping),
            n_mapped=len(targets),
            n_unmapped=len(unmapped),
            n_ambiguous=n_ambiguous,
            unmapped=unmapped[:MAX_NAMED_UNMAPPED],
        )
        logger.info(report.summary())
        return report
