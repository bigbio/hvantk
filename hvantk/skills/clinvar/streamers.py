"""ClinVar concrete streamers."""
from __future__ import annotations

from typing import Iterable

import hail as hl

from hvantk.core.streamers.variant_table import VariantTableStreamer


class ClinVarVariantTableStreamer(VariantTableStreamer):
    """ClinVar variant table streamer (info.CLNSIG, info.GENEINFO, info.CLNDN)."""

    def to_hail(self) -> hl.Table:
        return self._artifact.to_hail()

    def filter_to_genes(self, genes: Iterable[str]) -> hl.Table:
        ht = self.to_hail()
        gene_set = hl.literal(set(genes))
        gene_expr = hl.if_else(
            hl.is_defined(ht.info.GENEINFO) & (hl.len(ht.info.GENEINFO) > 0),
            ht.info.GENEINFO.split(":")[0],
            hl.missing(hl.tstr),
        )
        return ht.filter(gene_set.contains(gene_expr))

    def filter_by_pathogenicity(self, labels: Iterable[str]) -> hl.Table:
        ht = self.to_hail()
        label_set = hl.literal(set(labels))
        return ht.filter(ht.info.CLNSIG.any(lambda x: label_set.contains(x)))
