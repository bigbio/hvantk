"""PTM-variant analysis functions (Q1, Q3).

Implements:
- ptm_landscape: PTM-variant overlap counts and enrichment (Q1)
- ptm_population: Population-level allele frequency analysis at PTM sites (Q3)
- export_ptm_strata: Export PTM-stratified variant lists for external tools (e.g., PSROC)

Q2 (predictor evaluation) is composed at the workflow level: annotate variants
with PTM info (Phase 3), export strata, then run hvantk psroc independently.
"""

import json
import logging
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import hail as hl

from hvantk.core.constants import CLINVAR_PATHOGENIC_LABELS, CLINVAR_BENIGN_LABELS

logger = logging.getLogger(__name__)


@dataclass
class PTMLandscapeResult:
    """Result from PTM-variant landscape analysis (Q1)."""

    n_variants: int = 0
    n_pathogenic: int = 0
    n_benign: int = 0
    n_ptm_site_pathogenic: int = 0
    n_ptm_proximal_pathogenic: int = 0
    n_ptm_site_benign: int = 0
    n_ptm_proximal_benign: int = 0
    enrichment_odds_ratio: float = 0.0
    enrichment_ci_low: float = 0.0
    enrichment_ci_high: float = float("inf")
    enrichment_p_value: float = 1.0
    overlap_by_category: Dict[str, int] = field(default_factory=dict)
    category_enrichment: Dict[str, dict] = field(default_factory=dict)
    distance_distribution: Dict[int, int] = field(default_factory=dict)
    output_dir: str = ""

    def summary(self) -> str:
        lines = [
            "PTM-Variant Landscape:",
            f"  Variants: {self.n_variants:,} total, "
            f"{self.n_pathogenic:,} P/LP, {self.n_benign:,} B/LB",
            f"  P/LP at PTM site: {self.n_ptm_site_pathogenic:,}",
            f"  P/LP proximal: {self.n_ptm_proximal_pathogenic:,}",
            f"  Enrichment: OR={self.enrichment_odds_ratio:.2f}, "
            f"p={self.enrichment_p_value:.2e}",
        ]
        if self.overlap_by_category:
            lines.append("  By PTM category (P/LP):")
            for cat, n in sorted(self.overlap_by_category.items(), key=lambda x: -x[1]):
                lines.append(f"    {cat}: {n:,}")
        return "\n".join(lines)


@dataclass
class PTMPopulationResult:
    """Result from population-level PTM-variant analysis (Q3)."""

    n_variants: int = 0
    n_ptm_site: int = 0
    n_ptm_proximal: int = 0
    n_non_ptm: int = 0
    mean_af_ptm_site: float = 0.0
    mean_af_ptm_proximal: float = 0.0
    mean_af_non_ptm: float = 0.0
    n_zero_af_ptm: int = 0
    ptm_site_afs: List[float] = field(default_factory=list)
    proximal_afs: List[float] = field(default_factory=list)
    non_ptm_afs: List[float] = field(default_factory=list)
    ccr_mean_ptm: Optional[float] = None
    ccr_mean_non_ptm: Optional[float] = None
    output_dir: str = ""

    def summary(self) -> str:
        lines = [
            "PTM Population Analysis:",
            f"  Total variants: {self.n_variants:,}",
            f"  At PTM site: {self.n_ptm_site:,} (mean AF={self.mean_af_ptm_site:.2e})",
            f"  Proximal: {self.n_ptm_proximal:,} (mean AF={self.mean_af_ptm_proximal:.2e})",
            f"  Non-PTM: {self.n_non_ptm:,} (mean AF={self.mean_af_non_ptm:.2e})",
            f"  PTM sites with zero AF: {self.n_zero_af_ptm:,}",
        ]
        if self.ccr_mean_ptm is not None and self.ccr_mean_non_ptm is not None:
            lines.append(
                f"  Mean CCR: PTM={self.ccr_mean_ptm:.1f}, "
                f"non-PTM={self.ccr_mean_non_ptm:.1f}"
            )
        elif self.ccr_mean_ptm is not None:
            lines.append(f"  Mean CCR: PTM={self.ccr_mean_ptm:.1f}")
        elif self.ccr_mean_non_ptm is not None:
            lines.append(f"  Mean CCR: non-PTM={self.ccr_mean_non_ptm:.1f}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _extract_clnsig(ht: hl.Table) -> hl.Expression:
    """Extract scalar CLNSIG expression, handling array vs string dtype."""
    clnsig_dtype = ht.info.CLNSIG.dtype
    if isinstance(clnsig_dtype, hl.tarray):
        return hl.or_missing(
            hl.is_defined(ht.info.CLNSIG) & (hl.len(ht.info.CLNSIG) > 0),
            ht.info.CLNSIG[0],
        )
    return ht.info.CLNSIG


def _label_clinvar(ht: hl.Table) -> hl.Table:
    """Add _label field ('P', 'B', or 'other') based on ClinVar CLNSIG."""
    clnsig = _extract_clnsig(ht)
    p_labels = hl.literal(CLINVAR_PATHOGENIC_LABELS)
    b_labels = hl.literal(CLINVAR_BENIGN_LABELS)
    return ht.annotate(
        _label=hl.case()
        .when(hl.is_defined(clnsig) & p_labels.contains(clnsig), "P")
        .when(hl.is_defined(clnsig) & b_labels.contains(clnsig), "B")
        .default("other")
    )


def _odds_ratio_with_ci(table, confidence_level=0.95):
    """Compute OR, CI, and p-value from a 2x2 contingency table.

    Uses scipy.stats.contingency.odds_ratio (exact conditional) when
    available (scipy >= 1.8), otherwise falls back to Fisher's exact
    test for OR/p and Woolf's logit method for CI.

    Parameters
    ----------
    table : list of list
        2x2 contingency table [[a, b], [c, d]].
    confidence_level : float
        Confidence level for CI (default: 0.95).

    Returns
    -------
    tuple of (odds_ratio, ci_low, ci_high, p_value)
    """
    import math
    from scipy.stats import fisher_exact

    try:
        from scipy.stats.contingency import odds_ratio as scipy_or

        result = scipy_or(table, kind="conditional")
        ci = result.confidence_interval(confidence_level)
        _, p_value = fisher_exact(table)
        logger.debug("Using scipy.stats.contingency.odds_ratio for CI")
        return float(result.statistic), float(ci.low), float(ci.high), float(p_value)
    except ImportError:
        logger.debug("scipy odds_ratio not available; falling back to Fisher/Woolf")
    except Exception as exc:
        logger.debug(
            "scipy odds_ratio failed (%s); falling back to Fisher/Woolf",
            type(exc).__name__,
        )

    # Fallback: Fisher p-value + Woolf logit OR and CI (consistent estimator)
    _, p_value = fisher_exact(table)
    a, b = table[0]
    c, d = table[1]
    # Apply 0.5 continuity correction if any cell is zero
    if 0 in (a, b, c, d):
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    try:
        log_or = math.log(a * d / (b * c))
        or_val = math.exp(log_or)
        se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        from scipy.stats import norm

        z = norm.ppf(1 - (1 - confidence_level) / 2)
        ci_low = math.exp(log_or - z * se)
        ci_high = math.exp(log_or + z * se)
    except (ValueError, ZeroDivisionError):
        or_val = 0.0
        ci_low, ci_high = 0.0, float("inf")
    return float(or_val), ci_low, ci_high, float(p_value)


def export_ptm_strata(
    annotated_ht: hl.Table,
    output_dir: str,
) -> Dict[str, str]:
    """Export PTM-stratified variant lists from a PTM-annotated table.

    Splits a variant table (annotated by annotate_variants_with_ptm) into
    PTM-site and non-PTM variant lists in chr:pos:ref:alt format, suitable
    for downstream tools like ``hvantk psroc``.

    Parameters
    ----------
    annotated_ht : hl.Table
        Variant table with is_ptm_site and is_ptm_proximal fields
        (output of annotate_variants_with_ptm).
    output_dir : str
        Directory to write variant list files.

    Returns
    -------
    dict
        Maps stratum name to output file path:
        {"ptm": "<dir>/ptm_variants.txt", "non_ptm": "<dir>/non_ptm_variants.txt"}
    """
    os.makedirs(output_dir, exist_ok=True)

    strata = {
        "ptm": annotated_ht.filter(
            annotated_ht.is_ptm_site | annotated_ht.is_ptm_proximal
        ),
        "non_ptm": annotated_ht.filter(
            ~annotated_ht.is_ptm_site & ~annotated_ht.is_ptm_proximal
        ),
    }

    paths = {}
    for name, ht in strata.items():
        out_path = os.path.join(output_dir, f"{name}_variants.txt")
        ht_out = ht.annotate(
            _vid=hl.delimit(
                [
                    hl.str(ht.locus.contig),
                    hl.str(ht.locus.position),
                    ht.alleles[0],
                    ht.alleles[1],
                ],
                ":",
            )
        )
        ht_out.key_by().select("_vid").export(out_path, header=False)
        n = ht.count()
        logger.info(f"Exported {n:,} {name} variants to {out_path}")
        paths[name] = out_path

    return paths


# ---------------------------------------------------------------------------
# Q1: PTM-Variant Landscape
# ---------------------------------------------------------------------------


def ptm_landscape(
    clinvar_ht: hl.Table,
    ptm_ht: hl.Table,
    output_dir: str,
    flanking_codons: int = 5,
) -> PTMLandscapeResult:
    """PTM-variant overlap and enrichment analysis (Q1).

    Cross-references ClinVar P/LP variants with PTM sites. Computes overlap
    counts per PTM category, overall enrichment (Fisher's exact test), and
    distance distribution for plotting.

    Parameters
    ----------
    clinvar_ht : hl.Table
        ClinVar Hail Table keyed by (locus, alleles) with info.CLNSIG.
    ptm_ht : hl.Table
        PTM sites Hail Table from create_ptm_sites_tb.
    output_dir : str
        Directory for output files (landscape_summary.json).
    flanking_codons : int
        Flanking codons for proximal window (default: 5).

    Returns
    -------
    PTMLandscapeResult
    """
    from hvantk.algorithms.ptm.annotate import annotate_variants_with_ptm

    os.makedirs(output_dir, exist_ok=True)
    logger.info("Running PTM-variant landscape analysis (Q1)")

    # Annotate ClinVar with PTM info and assign P/B labels
    annotated = annotate_variants_with_ptm(clinvar_ht, ptm_ht, flanking_codons)
    annotated = _label_clinvar(annotated)

    # Aggregate counts in a single pass
    stats = annotated.aggregate(
        hl.struct(
            n_total=hl.agg.count(),
            n_path=hl.agg.filter(annotated._label == "P", hl.agg.count()),
            n_benign=hl.agg.filter(annotated._label == "B", hl.agg.count()),
            n_path_site=hl.agg.filter(
                (annotated._label == "P") & annotated.is_ptm_site, hl.agg.count()
            ),
            n_path_prox=hl.agg.filter(
                (annotated._label == "P") & annotated.is_ptm_proximal,
                hl.agg.count(),
            ),
            n_benign_site=hl.agg.filter(
                (annotated._label == "B") & annotated.is_ptm_site, hl.agg.count()
            ),
            n_benign_prox=hl.agg.filter(
                (annotated._label == "B") & annotated.is_ptm_proximal,
                hl.agg.count(),
            ),
        )
    )

    # Fisher's exact test with CI: (P/B) x (PTM/non-PTM)
    a = stats.n_path_site + stats.n_path_prox
    b = stats.n_benign_site + stats.n_benign_prox
    c = stats.n_path - a
    d = stats.n_benign - b
    odds_ratio, ci_low, ci_high, p_value = _odds_ratio_with_ci([[a, b], [c, d]])

    # Per-category overlap counts (P/LP at PTM sites, by category)
    ptm_path = annotated.filter(
        (annotated._label == "P") & (annotated.is_ptm_site | annotated.is_ptm_proximal)
    )
    ptm_path_exp = ptm_path.explode("ptm_types")
    category_counts = ptm_path_exp.aggregate(hl.agg.counter(ptm_path_exp.ptm_types))
    category_counts = {k: v for k, v in category_counts.items() if k is not None}

    # Per-category enrichment (Fisher's test per PTM type)
    # For each category: 2x2 table of (P/B) x (this_category / not_this_category)
    cat_benign = annotated.filter(
        (annotated._label == "B") & (annotated.is_ptm_site | annotated.is_ptm_proximal)
    )
    cat_benign_exp = cat_benign.explode("ptm_types")
    benign_cat_counts = cat_benign_exp.aggregate(
        hl.agg.counter(cat_benign_exp.ptm_types)
    )
    benign_cat_counts = {k: v for k, v in benign_cat_counts.items() if k is not None}

    category_enrichment = {}
    total_ptm_path = a  # P/LP at PTM (site + proximal)
    total_ptm_benign = b  # B/LB at PTM (site + proximal)
    for cat in sorted(set(category_counts) | set(benign_cat_counts)):
        n_path_cat = category_counts.get(cat, 0)
        n_benign_cat = benign_cat_counts.get(cat, 0)
        # 2x2: (P in cat, B in cat) vs (P not in cat, B not in cat)
        cat_or, cat_ci_lo, cat_ci_hi, cat_p = _odds_ratio_with_ci(
            [
                [n_path_cat, n_benign_cat],
                [
                    max(total_ptm_path - n_path_cat, 0),
                    max(total_ptm_benign - n_benign_cat, 0),
                ],
            ]
        )
        category_enrichment[cat] = {
            "n_pathogenic": n_path_cat,
            "n_benign": n_benign_cat,
            "odds_ratio": float(cat_or),
            "ci_low": float(cat_ci_lo),
            "ci_high": float(cat_ci_hi),
            "p_value": float(cat_p),
        }

    # Distance distribution (P/LP variants near PTM sites)
    near_ptm = annotated.filter(
        (annotated._label == "P") & hl.is_defined(annotated.ptm_distance)
    )
    distance_dist = near_ptm.aggregate(hl.agg.counter(near_ptm.ptm_distance))
    distance_dist = {int(k): v for k, v in distance_dist.items() if k is not None}

    result = PTMLandscapeResult(
        n_variants=stats.n_total,
        n_pathogenic=stats.n_path,
        n_benign=stats.n_benign,
        n_ptm_site_pathogenic=stats.n_path_site,
        n_ptm_proximal_pathogenic=stats.n_path_prox,
        n_ptm_site_benign=stats.n_benign_site,
        n_ptm_proximal_benign=stats.n_benign_prox,
        enrichment_odds_ratio=float(odds_ratio),
        enrichment_ci_low=float(ci_low),
        enrichment_ci_high=float(ci_high),
        enrichment_p_value=float(p_value),
        overlap_by_category=category_counts,
        category_enrichment=category_enrichment,
        distance_distribution=distance_dist,
        output_dir=output_dir,
    )

    # Write JSON summary
    with open(os.path.join(output_dir, "landscape_summary.json"), "w") as f:
        json.dump(
            {
                "n_variants": result.n_variants,
                "n_pathogenic": result.n_pathogenic,
                "n_benign": result.n_benign,
                "ptm_site": {
                    "pathogenic": result.n_ptm_site_pathogenic,
                    "benign": result.n_ptm_site_benign,
                },
                "ptm_proximal": {
                    "pathogenic": result.n_ptm_proximal_pathogenic,
                    "benign": result.n_ptm_proximal_benign,
                },
                "enrichment": {
                    "odds_ratio": result.enrichment_odds_ratio,
                    "ci_low": result.enrichment_ci_low,
                    "ci_high": result.enrichment_ci_high,
                    "p_value": result.enrichment_p_value,
                },
                "overlap_by_category": result.overlap_by_category,
                "category_enrichment": result.category_enrichment,
                "distance_distribution": {
                    str(k): v for k, v in result.distance_distribution.items()
                },
            },
            f,
            indent=2,
        )

    logger.info(result.summary())
    return result


# ---------------------------------------------------------------------------
# Q3: Population-Level PTM-Variant Burden
# ---------------------------------------------------------------------------


def ptm_population(
    gnomad_ht: hl.Table,
    ptm_ht: hl.Table,
    output_dir: str,
    ccr_ht: Optional[hl.Table] = None,
    af_field: str = "AF",
    flanking_codons: int = 5,
) -> PTMPopulationResult:
    """Population-level allele frequency analysis at PTM sites (Q3).

    Annotates gnomAD variants with PTM proximity and compares allele
    frequency distributions between PTM-site and non-PTM coding variants.
    This is where Hail scalability matters (~700M gnomAD variants).

    Parameters
    ----------
    gnomad_ht : hl.Table
        gnomAD variant table with allele frequency field.
    ptm_ht : hl.Table
        PTM sites Hail Table from create_ptm_sites_tb.
    output_dir : str
        Directory for output files (population_summary.json).
    ccr_ht : hl.Table, optional
        CCR table for constrained region cross-reference.
    af_field : str
        Name of the allele frequency field in gnomad_ht (default: "AF").
    flanking_codons : int
        Flanking codons for proximal window (default: 5).

    Returns
    -------
    PTMPopulationResult
    """
    from hvantk.algorithms.ptm.annotate import annotate_variants_with_ptm

    os.makedirs(output_dir, exist_ok=True)
    logger.info("Running PTM population analysis (Q3)")

    # Annotate gnomAD with PTM info
    annotated = annotate_variants_with_ptm(gnomad_ht, ptm_ht, flanking_codons)

    # Optionally annotate with CCR
    if ccr_ht is not None:
        annotated = annotated.annotate(_ccr_pct=ccr_ht[annotated.locus].ccr_pct)

    af = annotated[af_field]

    # Compute AF statistics in a single pass
    is_ptm = annotated.is_ptm_site
    is_prox = annotated.is_ptm_proximal
    is_non_ptm = ~is_ptm & ~is_prox

    stats = annotated.aggregate(
        hl.struct(
            n_total=hl.agg.count(),
            n_ptm_site=hl.agg.filter(is_ptm, hl.agg.count()),
            n_ptm_prox=hl.agg.filter(is_prox, hl.agg.count()),
            n_non_ptm=hl.agg.filter(is_non_ptm, hl.agg.count()),
            af_ptm_site=hl.agg.filter(is_ptm, hl.agg.stats(af)),
            af_ptm_prox=hl.agg.filter(is_prox, hl.agg.stats(af)),
            af_non_ptm=hl.agg.filter(is_non_ptm, hl.agg.stats(af)),
            n_zero_af_ptm=hl.agg.filter(is_ptm & (af == 0), hl.agg.count()),
            afs_ptm_site=hl.agg.filter(is_ptm, hl.agg.take(af, 5000)),
            afs_ptm_prox=hl.agg.filter(is_prox, hl.agg.take(af, 5000)),
            afs_non_ptm=hl.agg.filter(is_non_ptm, hl.agg.take(af, 5000)),
        )
    )

    # CCR comparison
    ccr_ptm = None
    ccr_non_ptm = None
    if ccr_ht is not None:
        ccr_stats = annotated.aggregate(
            hl.struct(
                ptm=hl.agg.filter(
                    is_ptm & hl.is_defined(annotated._ccr_pct),
                    hl.agg.stats(annotated._ccr_pct),
                ),
                non_ptm=hl.agg.filter(
                    is_non_ptm & hl.is_defined(annotated._ccr_pct),
                    hl.agg.stats(annotated._ccr_pct),
                ),
            )
        )
        if ccr_stats.ptm.n > 0:
            ccr_ptm = ccr_stats.ptm.mean
        if ccr_stats.non_ptm.n > 0:
            ccr_non_ptm = ccr_stats.non_ptm.mean

    result = PTMPopulationResult(
        n_variants=stats.n_total,
        n_ptm_site=stats.n_ptm_site,
        n_ptm_proximal=stats.n_ptm_prox,
        n_non_ptm=stats.n_non_ptm,
        mean_af_ptm_site=(stats.af_ptm_site.mean if stats.af_ptm_site.n > 0 else 0.0),
        mean_af_ptm_proximal=(
            stats.af_ptm_prox.mean if stats.af_ptm_prox.n > 0 else 0.0
        ),
        mean_af_non_ptm=(stats.af_non_ptm.mean if stats.af_non_ptm.n > 0 else 0.0),
        n_zero_af_ptm=stats.n_zero_af_ptm,
        ptm_site_afs=[float(x) for x in stats.afs_ptm_site],
        proximal_afs=[float(x) for x in stats.afs_ptm_prox],
        non_ptm_afs=[float(x) for x in stats.afs_non_ptm],
        ccr_mean_ptm=ccr_ptm,
        ccr_mean_non_ptm=ccr_non_ptm,
        output_dir=output_dir,
    )

    # Write JSON summary
    with open(os.path.join(output_dir, "population_summary.json"), "w") as f:
        json.dump(
            {
                "n_variants": result.n_variants,
                "n_ptm_site": result.n_ptm_site,
                "n_ptm_proximal": result.n_ptm_proximal,
                "n_non_ptm": result.n_non_ptm,
                "mean_af": {
                    "ptm_site": result.mean_af_ptm_site,
                    "ptm_proximal": result.mean_af_ptm_proximal,
                    "non_ptm": result.mean_af_non_ptm,
                },
                "n_zero_af_ptm": result.n_zero_af_ptm,
                "ptm_site_afs": result.ptm_site_afs,
                "proximal_afs": result.proximal_afs,
                "non_ptm_afs": result.non_ptm_afs,
                "ccr_mean_ptm": result.ccr_mean_ptm,
                "ccr_mean_non_ptm": result.ccr_mean_non_ptm,
            },
            f,
            indent=2,
        )

    logger.info(result.summary())
    return result
