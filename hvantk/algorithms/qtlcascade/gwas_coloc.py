"""GWAS → effector colocalization (FinnGen GWAS × eQTL Catalogue cis-eQTL).

ABF colocalization that ranks the cis genes at a GWAS locus by the posterior
probability that the GWAS signal and the gene's cis-eQTL share a causal variant
(H4). All summary statistics are streamed via **remote tabix** — no bulk
downloads. The single-variant ABF here is the fast screen; for rigour it should
be confirmed by fine-mapping (see :mod:`hvantk.algorithms.qtlcascade.finemap`),
which separates genuine colocalization from single-variant-ABF artifacts.

Reuses the validated Wakefield kernel (:func:`coloc.compute_log_abf`) with
trait-specific prior variances (case-control GWAS vs quantitative eQTL).

References
----------
- Giambartolomei et al. (2014) PLoS Genet — coloc
- Wakefield (2009) Am J Hum Genet — ABF
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import pysam

from hvantk.algorithms.qtlcascade.coloc import (
    compute_log_abf,
    _logsumexp,
    _logdiff,
)
from hvantk.algorithms.qtlcascade.constants import (
    DEFAULT_COLOC_P1,
    DEFAULT_COLOC_P2,
    DEFAULT_COLOC_P12,
    DEFAULT_GWAS_W_CC,
    DEFAULT_EQTL_W_QUANT,
    DEFAULT_COLOC_MIN_SNPS,
    DEFAULT_COLOC_WINDOW_KB,
    FINNGEN_R10_SUMSTATS_URL,
    EQTL_CATALOGUE_SUMSTATS_URL,
    EQTL_CATALOGUE_DEFAULT_STUDY,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Remote-tabix helpers
# ---------------------------------------------------------------------------


def _ensure_ca_bundle() -> None:
    """htslib remote-TLS needs a CA bundle; point it at certifi if unset."""
    if not os.environ.get("CURL_CA_BUNDLE"):
        try:
            import certifi

            os.environ["CURL_CA_BUNDLE"] = certifi.where()
        except Exception as exc:  # pragma: no cover - certifi is a hard dep
            logger.debug("certifi unavailable for CA bundle: %s", exc)


def _contig(tb: "pysam.TabixFile", chrom: str) -> str:
    """Resolve a contig name against a tabix file's naming (with/without 'chr')."""
    if chrom in tb.contigs:
        return chrom
    alt = f"chr{chrom}" if not chrom.startswith("chr") else chrom[3:]
    return alt if alt in tb.contigs else chrom


def finngen_url(endpoint: str) -> str:
    return FINNGEN_R10_SUMSTATS_URL.format(endpoint=endpoint)


def eqtl_catalogue_url(dataset: str, study: str = EQTL_CATALOGUE_DEFAULT_STUDY) -> str:
    """Build an eQTL Catalogue URL from a dataset id (e.g. ``QTD000251``).

    If ``dataset`` already looks like a URL/path it is returned unchanged.
    """
    if "://" in dataset or dataset.endswith(".tsv.gz") or os.path.exists(dataset):
        return dataset
    return EQTL_CATALOGUE_SUMSTATS_URL.format(study=study, dataset=dataset)


def fetch_finngen_region(
    endpoint: str, chrom: str, start: int, end: int
) -> dict[tuple[int, str, str], tuple[float, float, float]]:
    """FinnGen R10 GWAS over a region → ``{(pos, ref, alt): (beta, se, pval)}``.

    Effects are ALT-referenced. Cols: chrom pos ref alt rsids nearest_genes
    pval mlogp beta sebeta af*.
    """
    _ensure_ca_bundle()
    out: dict[tuple[int, str, str], tuple[float, float, float]] = {}
    with pysam.TabixFile(finngen_url(endpoint)) as tb:
        for rec in tb.fetch(_contig(tb, chrom), start, end):
            f = rec.split("\t")
            try:
                se = float(f[9])
                if se > 0:
                    out[(int(f[1]), f[2], f[3])] = (float(f[8]), se, float(f[6]))
            except (ValueError, IndexError):
                continue
    return out


def fetch_eqtl_region(
    dataset: str, chrom: str, start: int, end: int, study: str = EQTL_CATALOGUE_DEFAULT_STUDY
) -> dict[str, list[tuple[tuple[int, str, str], float, float, float]]]:
    """eQTL Catalogue cis-eQTL over a region → ``{gene_id: [(key, beta, se, pval), ...]}``.

    Effects are ALT-referenced. Cols: gene_id chrom pos ref alt variant
    ma_samples maf pvalue beta se. ``dataset`` may be a QTD id, a full URL, or a
    local tabix path.
    """
    _ensure_ca_bundle()
    out: dict[str, list[tuple[tuple[int, str, str], float, float, float]]] = {}
    with pysam.TabixFile(eqtl_catalogue_url(dataset, study)) as tb:
        for rec in tb.fetch(_contig(tb, chrom), start, end):
            f = rec.split("\t")
            try:
                se = float(f[10])
                if se <= 0:  # guard against div-by-zero in the ABF kernel
                    continue
                key = (int(f[2]), f[3], f[4])
                out.setdefault(f[0].split(".")[0], []).append(
                    (key, float(f[9]), se, float(f[8]))
                )
            except (ValueError, IndexError):
                continue
    return out


# ---------------------------------------------------------------------------
# ABF coloc (two traits, trait-specific W) — reuses the validated kernel
# ---------------------------------------------------------------------------


def coloc_abf_two_traits(
    beta1: np.ndarray,
    se1: np.ndarray,
    beta2: np.ndarray,
    se2: np.ndarray,
    W1: float = DEFAULT_GWAS_W_CC,
    W2: float = DEFAULT_EQTL_W_QUANT,
    p1: float = DEFAULT_COLOC_P1,
    p2: float = DEFAULT_COLOC_P2,
    p12: float = DEFAULT_COLOC_P12,
) -> dict:
    """Wakefield/Giambartolomei ABF coloc with **trait-specific** prior variances.

    Identical hypothesis algebra to :func:`coloc.coloc_abf`, but uses W1 for
    trait 1 (e.g. case-control GWAS) and W2 for trait 2 (e.g. quantitative
    eQTL). Returns H0–H4 posteriors + ``n_variants`` + ``lead_idx``.
    """
    n = len(beta1)
    if n == 0:
        return {"H0": 1.0, "H1": 0.0, "H2": 0.0, "H3": 0.0, "H4": 0.0,
                "n_variants": 0, "lead_idx": -1}
    la1 = compute_log_abf(np.asarray(beta1), np.asarray(se1), W1)
    la2 = compute_log_abf(np.asarray(beta2), np.asarray(se2), W2)
    s1, s2 = _logsumexp(la1), _logsumexp(la2)
    la_both = la1 + la2
    s_both = _logsumexp(la_both)
    log_h = np.array([
        0.0,
        np.log(p1) + s1,
        np.log(p2) + s2,
        np.log(p1) + np.log(p2) + _logdiff(s1 + s2, s_both),
        np.log(p12) + s_both,
    ])
    post = np.exp(log_h - np.max(log_h))
    post /= post.sum()
    return {"H0": float(post[0]), "H1": float(post[1]), "H2": float(post[2]),
            "H3": float(post[3]), "H4": float(post[4]),
            "n_variants": n, "lead_idx": int(np.argmax(la_both))}


# ---------------------------------------------------------------------------
# Region driver
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GwasColocResult:
    """Ranked coloc result for one GWAS locus × one eQTL dataset."""

    table: pd.DataFrame          # gene_id, n_snp, eqtl_minp, PP3, PP4 (PP4-sorted)
    gwas_min_p: float
    n_genes_tested: int
    region: str


def run_locus_coloc(
    *,
    gwas: dict[tuple[int, str, str], tuple[float, float, float]],
    eqtl: dict[str, list[tuple[tuple[int, str, str], float, float, float]]],
    region: str,
    min_snps: int = DEFAULT_COLOC_MIN_SNPS,
    W1: float = DEFAULT_GWAS_W_CC,
    W2: float = DEFAULT_EQTL_W_QUANT,
    p1: float = DEFAULT_COLOC_P1,
    p2: float = DEFAULT_COLOC_P2,
    p12: float = DEFAULT_COLOC_P12,
    gene_symbols: Optional[dict[str, str]] = None,
) -> GwasColocResult:
    """ABF coloc of a GWAS region against every overlapping cis gene.

    Variants are matched on position and allele *set* — ref/alt may be stored in
    either order across sources (the ABF kernel uses z^2, so allele orientation
    does not affect H4). Pure in-memory — fetch with :func:`fetch_finngen_region`
    / :func:`fetch_eqtl_region` first.
    """
    gwas_min_p = min((v[2] for v in gwas.values()), default=float("nan"))
    rows = []
    for ensg, recs in eqtl.items():
        b1, s1, b2, s2, eqtl_ps = [], [], [], [], []
        for (pos, eref, ealt), ebeta, ese, epval in recs:
            # Match order-independently: sources may store ref/alt swapped.
            # The ABF kernel uses z^2 (sign-independent), so no beta flip needed.
            g = gwas.get((pos, eref, ealt)) or gwas.get((pos, ealt, eref))
            if g is None:
                continue
            b1.append(g[0]); s1.append(g[1])
            b2.append(ebeta); s2.append(ese); eqtl_ps.append(epval)
        if len(b1) < min_snps:
            continue
        res = coloc_abf_two_traits(
            np.array(b1), np.array(s1), np.array(b2), np.array(s2),
            W1=W1, W2=W2, p1=p1, p2=p2, p12=p12,
        )
        rows.append({
            "gene_id": ensg,
            "gene": (gene_symbols or {}).get(ensg, ensg),
            "n_snp": len(b1),
            "eqtl_min_p": min(eqtl_ps),
            "PP3": res["H3"],
            "PP4": res["H4"],
        })
    cols = ["gene_id", "gene", "n_snp", "eqtl_min_p", "PP3", "PP4"]
    df = (pd.DataFrame(rows, columns=cols).sort_values("PP4", ascending=False)
          .reset_index(drop=True) if rows else pd.DataFrame(columns=cols))
    return GwasColocResult(table=df, gwas_min_p=gwas_min_p,
                           n_genes_tested=len(df), region=region)


def run_finngen_eqtl_coloc(
    *,
    endpoint: str,
    chrom: str,
    lead: int,
    eqtl_dataset: str,
    window_kb: int = DEFAULT_COLOC_WINDOW_KB,
    eqtl_study: str = EQTL_CATALOGUE_DEFAULT_STUDY,
    min_snps: int = DEFAULT_COLOC_MIN_SNPS,
    W1: float = DEFAULT_GWAS_W_CC,
    W2: float = DEFAULT_EQTL_W_QUANT,
    p1: float = DEFAULT_COLOC_P1,
    p2: float = DEFAULT_COLOC_P2,
    p12: float = DEFAULT_COLOC_P12,
    gene_symbols: Optional[dict[str, str]] = None,
) -> GwasColocResult:
    """End-to-end ABF coloc: fetch FinnGen × eQTL Catalogue region, rank effectors."""
    half = window_kb * 1000
    start, end = max(0, lead - half), lead + half  # clamp left edge near contig start
    region = f"chr{chrom}:{start}-{end}"
    logger.info("coloc %s × %s @ %s", endpoint, eqtl_dataset, region)
    gwas = fetch_finngen_region(endpoint, chrom, start, end)
    eqtl = fetch_eqtl_region(eqtl_dataset, chrom, start, end, study=eqtl_study)
    logger.info("  GWAS variants=%d; eQTL genes=%d", len(gwas), len(eqtl))
    return run_locus_coloc(
        gwas=gwas, eqtl=eqtl, region=region, min_snps=min_snps,
        W1=W1, W2=W2, p1=p1, p2=p2, p12=p12, gene_symbols=gene_symbols,
    )
