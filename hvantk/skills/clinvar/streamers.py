"""ClinVar concrete streamers."""
from __future__ import annotations

from typing import Iterable, Optional, Sequence, Set

import hail as hl

from hvantk.core.streamers.variant_table import VariantTableStreamer
from hvantk.skills.clinvar.shared.constants import (
    CLINVAR_PATHOGENIC_LABELS,
    CLINVAR_BENIGN_LABELS,
)


def apply_clinvar_training_labels(
    ht: hl.Table,
    *,
    pathogenic_labels: Sequence[str],
    benign_labels: Sequence[str],
    gene_set: Optional[Iterable[str]] = None,
    disease_terms: Optional[Iterable[str]] = None,
    label_column: str = "rf_label",
) -> hl.Table:
    """Add a TP/TN training-label column derived from ClinVar CLNSIG/CLNDN.

    Ports the TP/TN rule from ``ClinvarDataStreamer``
    (skills/clinvar/pipelines/training_set.py) to a whole-table function:

      - derive a ``gene`` column from ``info.GENEINFO`` (first token before
        ``':'``) if absent;
      - exclude synonymous variants (``info.MC`` consequence
        ``'synonymous_variant'``), matching the legacy ``stream()`` pre-filter;
      - gene-based TP: ``info.CLNSIG`` any in ``pathogenic_labels`` AND
        (``gene_set`` empty OR ``gene`` in ``gene_set``);
      - disease-based TP: any ``CLNDN`` token (split on ``'|'``, spaces->``'_'``,
        lowercased) in ``disease_terms``;
      - TN: ``info.CLNSIG`` any in ``benign_labels``;
      - keep rows where (TP) xor (TN); set ``label_column`` to ``'TP'``/``'TN'``
        and drop the rest.

    Preserves the variant keys (locus, alleles) and the ``gene`` column for
    downstream gene-based annotation.

    Note: this assumes ``info.CLNSIG`` is an array field (it calls
    ``.any(...)`` on it), matching the legacy ``ClinvarDataStreamer``. This
    intentionally differs from
    ``ClinVarVariantTableStreamer.filter_by_pathogenicity``, which handles both
    array and scalar CLNSIG.
    """
    gene_set = set(gene_set) if gene_set else set()
    # Normalize disease terms the same way the legacy streamer did.
    normalized_disease_terms: Set[str] = {
        d.replace(" ", "_").lower() for d in (disease_terms or set())
    }

    info_fields = ht.row.dtype["info"].fields

    # --- derive `gene` from info.GENEINFO if absent (legacy setup()) ---
    if "gene" not in ht.row:
        if "GENEINFO" in info_fields:
            ht = ht.annotate(
                gene=hl.if_else(
                    hl.is_defined(ht.info.GENEINFO)
                    & (hl.len(ht.info.GENEINFO) > 0),
                    ht.info.GENEINFO.split(":")[0],
                    hl.missing(hl.tstr),
                )
            )
        else:
            ht = ht.annotate(gene=hl.missing(hl.tstr))

    # --- consequence array for the synonymous pre-filter (legacy setup()) ---
    if "MC" in info_fields:
        consequence = ht.info.MC.map(
            lambda x: hl.if_else(
                hl.is_defined(x) & (hl.len(x.split("|")) > 1),
                x.split("|")[1],
                "",
            )
        )
    else:
        consequence = hl.empty_array(hl.tstr)
    ht = ht.annotate(_consequence=consequence)

    # Exclude synonymous variants (legacy stream() pre-filter).
    ht = ht.filter(
        hl.is_missing(ht._consequence)
        | ~ht._consequence.any(lambda x: x == "synonymous_variant")
    )

    # --- disease-based TP (legacy process_chunk) ---
    if normalized_disease_terms:
        disease_set = hl.literal(normalized_disease_terms)
        clndn = hl.or_else(ht.info.CLNDN, "")
        tokens = (
            hl.str(clndn).split(r"\|").map(lambda t: t.replace(" ", "_").lower())
        )
        disease_tp = tokens.any(lambda t: disease_set.contains(t))
    else:
        disease_tp = hl.literal(False)

    # --- gene-based TP (legacy process_chunk) ---
    gene_filter = (
        hl.literal(True)
        if not gene_set
        else hl.literal(gene_set).contains(ht.gene)
    )
    gene_tp = (
        ht.info.CLNSIG.any(lambda x: hl.set(pathogenic_labels).contains(x))
        & gene_filter
    )

    is_tp_site = gene_tp | disease_tp
    is_tn_site = ht.info.CLNSIG.any(lambda x: hl.set(benign_labels).contains(x))

    ht = ht.annotate(_is_tp_site=is_tp_site, _is_tn_site=is_tn_site)
    # Keep rows where exactly one of TP/TN holds.
    ht = ht.filter(ht._is_tp_site != ht._is_tn_site)
    ht = ht.annotate(
        **{
            label_column: hl.case()
            .when(ht._is_tp_site, "TP")
            .when(ht._is_tn_site, "TN")
            .or_missing()
        }
    )
    ht = ht.filter(hl.is_defined(ht[label_column]))
    return ht.drop("_consequence", "_is_tp_site", "_is_tn_site")


class ClinVarVariantTableStreamer(VariantTableStreamer):
    """ClinVar variant table streamer (info.CLNSIG, info.GENEINFO, info.CLNDN)."""

    PATHOGENIC_LABELS = CLINVAR_PATHOGENIC_LABELS
    BENIGN_LABELS = CLINVAR_BENIGN_LABELS

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
        clnsig = ht.info.CLNSIG
        # CLNSIG may be imported as an array or a scalar string depending on
        # the VCF header's Number= field; handle both (mirrors the dtype
        # check in algorithms/ptm/analysis.py:_extract_clnsig).
        if isinstance(clnsig.dtype, hl.tarray):
            return ht.filter(
                hl.is_defined(clnsig)
                & clnsig.any(lambda x: label_set.contains(x))
            )
        return ht.filter(hl.is_defined(clnsig) & label_set.contains(clnsig))

    def label_training_set(
        self,
        *,
        gene_set: Optional[Iterable[str]] = None,
        disease_terms: Optional[Iterable[str]] = None,
        label_column: str = "rf_label",
    ) -> hl.Table:
        """Derive a TP/TN training-label column from this ClinVar table."""
        return apply_clinvar_training_labels(
            self.to_hail(),
            pathogenic_labels=self.PATHOGENIC_LABELS,
            benign_labels=self.BENIGN_LABELS,
            gene_set=gene_set,
            disease_terms=disease_terms,
            label_column=label_column,
        )
