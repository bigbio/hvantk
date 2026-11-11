"""ClinVar concrete streamers."""

from __future__ import annotations

from typing import Iterable, Sequence

import hail as hl

from hvantk.core.streamers.variant_table import VariantTableStreamer
from hvantk.skills.clinvar.shared.constants import (
    CLINVAR_BENIGN_LABELS,
    CLINVAR_PATHOGENIC_LABELS,
    CLINVAR_REVIEW_STARS,
)


def _joined(expr) -> hl.StringExpression:
    """A VCF INFO string, re-joined with ``,`` when ``import_vcf`` split it into an array.

    Fields declared ``Number=.`` (CLNSIG, CLNREVSTAT, CLNDN, ...) arrive as ``array<str>``
    split at every comma, including commas inside one value
    (``criteria_provided,_single_submitter``). An empty array or string is missing.
    """
    if isinstance(expr.dtype, hl.tarray):
        return hl.or_missing(hl.is_defined(expr) & (hl.len(expr) > 0), hl.delimit(expr, ","))
    return hl.or_missing(hl.is_defined(expr) & (expr != ""), expr)


def label_clinvar_variants(
    ht: hl.Table,
    *,
    pathogenic: Sequence[str],
    benign: Sequence[str],
    consequence: str = "missense_variant",
    snv_only: bool = True,
) -> hl.Table:
    """Label ClinVar variants 1 (pathogenic) / 0 (benign) with their review stars.

    Keeps the rows that pass the variant filters -- biallelic SNVs when ``snv_only``,
    and an ``MC`` that lists ``consequence`` -- and returns them keyed ``(locus,
    alleles)`` with:

    - ``vid``: ``chrom:pos:ref:alt``;
    - ``clinvar_clnsig``: CLNSIG as released (re-joined with ``,``);
    - ``clinvar_review_status``: CLNREVSTAT (re-joined) and ``clinvar_review_stars``
      (int32, 0-4, :data:`CLINVAR_REVIEW_STARS`; anything unlisted is 0);
    - ``consequence``: the distinct MC terms, sorted, joined with ``,`` (a malformed
      element without a ``|`` contributes nothing);
    - ``label``: int32 1 when the primary CLNSIG (text before the first ``|``) is in
      ``pathogenic``, 0 when it is in ``benign``, missing otherwise -- kept, not
      filtered, so a caller can count every step;
    - ``geneinfo_symbols`` / ``geneinfo_entrez``: the GENEINFO genes in listed order
      (an Entrez id is missing when a gene has none).

    A missing or empty CLNSIG, CLNREVSTAT, MC or GENEINFO gives missing values or an
    empty list, never an error.
    """
    mc = ht.info.MC
    if not isinstance(mc.dtype, hl.tarray):
        mc = hl.or_missing(hl.is_defined(mc), mc.split(","))
    # Each MC element is "SO:0001583|missense_variant"; Hail raises on an out-of-bounds
    # index rather than returning missing, so every split is guarded.
    mc_terms = hl.or_else(mc, hl.empty_array(hl.tstr)).map(
        lambda x: hl.if_else(
            hl.is_defined(x) & (hl.len(x.split(r"\|")) > 1), x.split(r"\|")[1], ""
        )
    )

    geneinfo = hl.or_else(_joined(ht.info.GENEINFO), "")
    ht = ht.annotate(
        _mc=mc_terms,
        _genes=hl.if_else(
            geneinfo == "", hl.empty_array(hl.tstr), geneinfo.split(r"\|")
        ),
    )

    keep = ht._mc.contains(consequence)
    if snv_only:
        keep = keep & hl.if_else(
            hl.len(ht.alleles) == 2, hl.is_snp(ht.alleles[0], ht.alleles[1]), False
        )
    ht = ht.filter(keep)

    clnsig = _joined(ht.info.CLNSIG)
    revstat = _joined(ht.info.CLNREVSTAT)
    primary = clnsig.split(r"\|")[0]
    ht = ht.annotate(
        vid=hl.format("%s:%d:", ht.locus.contig, ht.locus.position)
        + hl.delimit(ht.alleles, ":"),
        clinvar_clnsig=clnsig,
        clinvar_review_status=revstat,
        clinvar_review_stars=hl.or_else(
            hl.literal(CLINVAR_REVIEW_STARS, hl.tdict(hl.tstr, hl.tint32)).get(revstat),
            0,
        ),
        consequence=hl.delimit(
            hl.sorted(hl.array(hl.set(ht._mc.filter(lambda t: t != "")))), ","
        ),
        label=hl.case()
        .when(hl.literal(set(pathogenic), hl.tset(hl.tstr)).contains(primary), 1)
        .when(hl.literal(set(benign), hl.tset(hl.tstr)).contains(primary), 0)
        .or_missing(),
        geneinfo_symbols=ht._genes.map(lambda g: g.split(":")[0]),
        geneinfo_entrez=ht._genes.map(
            lambda g: hl.or_missing(hl.len(g.split(":")) > 1, g.split(":")[1])
        ),
    )
    return ht.select(
        "vid",
        "clinvar_clnsig",
        "clinvar_review_status",
        "clinvar_review_stars",
        "consequence",
        "label",
        "geneinfo_symbols",
        "geneinfo_entrez",
    )


class ClinVarVariantTableStreamer(VariantTableStreamer):
    """ClinVar variant table streamer (info.CLNSIG, info.CLNREVSTAT, info.MC, info.GENEINFO)."""

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
                hl.is_defined(clnsig) & clnsig.any(lambda x: label_set.contains(x))
            )
        return ht.filter(hl.is_defined(clnsig) & label_set.contains(clnsig))

    def label_training_set(self) -> hl.Table:
        """Label this ClinVar table 1/0 (see :func:`label_clinvar_variants`).

        Missense SNVs whose primary CLNSIG is pathogenic or benign, with their review
        stars; rows with neither label are dropped.
        """
        ht = label_clinvar_variants(
            self.to_hail(),
            pathogenic=self.PATHOGENIC_LABELS,
            benign=self.BENIGN_LABELS,
        )
        return ht.filter(hl.is_defined(ht.label))
