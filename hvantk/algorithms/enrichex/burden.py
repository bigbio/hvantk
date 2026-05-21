"""
Burden testing module using Hail-native regression.

This module provides functions for computing per-sample burden scores
for gene sets and testing their association with phenotypes using
Hail's native logistic and linear regression functions.

The implementation follows the pattern from alphagenome/logreg_burden_test.py,
which uses a two-step aggregation:
1. Variants → Genes (per sample)
2. Genes → Gene Sets (per sample)

This keeps all computation in Hail's distributed framework.
"""

from __future__ import annotations

import logging
import time
import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, List, Optional

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd

try:  # Optional dependency for burden analysis
    import hail as hl
except ModuleNotFoundError as exc:  # pragma: no cover - depends on env
    hl = None  # type: ignore
    _HAIL_IMPORT_ERROR = exc
else:
    _HAIL_IMPORT_ERROR = None

from hvantk.algorithms.enrichex.constants import (
    VARIANT_CLASS_PRESETS,
    _DEPRECATED_AGGREGATION_ALIASES,
)
from hvantk.core.models.backends import algorithm, Backend
from hvantk.core.utils.table_utils import field_exists, resolve_field

logger = logging.getLogger(__name__)


def _resolve_genotype_aggregation(method: str) -> str:
    """Resolve deprecated genotype aggregation aliases to current names.

    Parameters
    ----------
    method : str
        Genotype aggregation method name (may be deprecated alias).

    Returns
    -------
    str
        Resolved method name.
    """
    if method in _DEPRECATED_AGGREGATION_ALIASES:
        new_name = _DEPRECATED_AGGREGATION_ALIASES[method]
        warnings.warn(
            f"Genotype aggregation '{method}' is deprecated, use '{new_name}' instead. "
            f"'{method}' approximates compound heterozygosity by counting genes with "
            f">= 2 heterozygous qualifying variants; this over-counts when variants "
            f"are in cis. For true compound-het calling, use phased genotypes.",
            DeprecationWarning,
            stacklevel=3,
        )
        return new_name
    return method


def _require_hail() -> None:
    if hl is None:  # pragma: no cover - depends on env
        raise ImportError(
            "Hail is required for EnrichEx burden analysis. "
            "Install hvantk with the 'hail' extra."
        ) from _HAIL_IMPORT_ERROR


@dataclass
class VariantFilter:
    """Criteria for qualifying variants in burden analysis.

    Attributes
    ----------
    max_af : float
        Maximum allele frequency threshold
    min_score : Optional[float]
        Minimum prediction score threshold (e.g., CADD, REVEL).
        None disables the filter.
    consequences : Optional[List[str]]
        List of qualifying VEP consequences
    pass_only : bool
        Only include PASS variants
    min_gq : int
        Minimum genotype quality
    min_dp : int
        Minimum depth
    af_field : str
        Field name for allele frequency
    score_field : str
        Field name for prediction score (e.g., ``"cadd_phred"``,
        ``"vep.CADD_PHRED"``, ``"REVEL"``).  Supports dot notation
        for nested structs.
    consequence_field : str
        Field name for consequence annotation
    """

    max_af: float = 0.01
    min_score: Optional[float] = 20.0
    consequences: Optional[List[str]] = None
    pass_only: bool = True
    min_gq: int = 20
    min_dp: int = 10

    # Field name mappings (configurable for different schemas)
    af_field: str = "gnomad_af"
    score_field: str = "cadd_phred"
    consequence_field: str = "consequence"

    def to_hail_expr(
        self, mt: hl.MatrixTable, row_only: bool = True
    ) -> hl.expr.BooleanExpression:
        """Convert filter criteria to Hail boolean expression.

        Parameters
        ----------
        mt : hl.MatrixTable
            MatrixTable to filter
        row_only : bool
            If True (default), only include row-level filters (AF, score,
            consequence, PASS).  Entry-level filters (GQ, DP) are excluded
            because they cannot be used in ``filter_rows()``.  Set to False
            to include entry-level filters (for use with ``filter_entries()``).

        Returns
        -------
        hl.expr.BooleanExpression
            Combined filter expression
        """
        _require_hail()
        filters = []

        # Row-level filters
        # Allele frequency filter
        if self.max_af is not None and self.max_af < 1.0:
            if field_exists(mt.row, self.af_field):
                filters.append(resolve_field(mt, self.af_field) <= self.max_af)
            else:
                logger.warning(f"AF field '{self.af_field}' not found in MT")

        # Prediction score filter (CADD, REVEL, etc.)
        if self.min_score is not None:
            if field_exists(mt.row, self.score_field):
                filters.append(resolve_field(mt, self.score_field) >= self.min_score)
            else:
                logger.warning(f"Score field '{self.score_field}' not found in MT")

        # Consequence filter
        if self.consequences:
            if field_exists(mt.row, self.consequence_field):
                filters.append(
                    hl.literal(self.consequences).contains(
                        resolve_field(mt, self.consequence_field)
                    )
                )
            else:
                logger.warning(
                    f"Consequence field '{self.consequence_field}' not found in MT"
                )

        # PASS filter
        if self.pass_only:
            if "filters" in mt.row:
                filters.append(hl.len(mt.filters) == 0)

        # Entry-level filters (only when explicitly requested)
        if not row_only:
            if self.min_gq is not None and self.min_gq > 0 and "GQ" in mt.entry:
                filters.append(mt.GQ >= self.min_gq)

            if self.min_dp is not None and self.min_dp > 0 and "DP" in mt.entry:
                filters.append(mt.DP >= self.min_dp)

        # Combine all filters
        if filters:
            return hl.all(*filters)
        else:
            return hl.bool(True)


def compute_per_gene_burden_mt(
    mt: "hl.MatrixTable",
    gene_field: str = "SYMBOL",
    variant_filter: Optional["VariantFilter"] = None,
) -> "hl.MatrixTable":
    """Aggregate variants to per-gene burden per sample.

    This is the first step of burden computation (variants -> genes).
    The returned MatrixTable has rows keyed by gene, columns keyed by
    sample, and three entry fields: ``hets`` (int), ``homs`` (int),
    ``multi_het`` (bool).

    Parameters
    ----------
    mt : hl.MatrixTable
        Annotated variant MatrixTable.  Expected to be pre-filtered.
    gene_field : str
        Row field containing gene symbol.
    variant_filter : VariantFilter, optional
        Optional variant-class filter (e.g., restrict to LoF).

    Returns
    -------
    hl.MatrixTable
        Per-gene burden MT (rows=genes, cols=samples).
    """
    _require_hail()

    if gene_field not in mt.row:
        raise ValueError(f"Gene field '{gene_field}' not found in MatrixTable")

    n_variants_in, n_samples_in = mt.count()
    logger.info(
        "compute_per_gene_burden_mt: input %d variants, %d samples",
        n_variants_in,
        n_samples_in,
    )

    # Apply variant filter if provided
    if variant_filter is not None:
        filter_expr = variant_filter.to_hail_expr(mt)
        mt = mt.filter_rows(filter_expr)
        n_after_filter = mt.count_rows()
        logger.info(
            "  %d variants after variant_filter (%d removed)",
            n_after_filter,
            n_variants_in - n_after_filter,
        )

    # Filter to rows with gene annotation
    mt = mt.filter_rows(hl.is_defined(mt[gene_field]))

    # Aggregate variants -> genes per sample
    mt_genes = mt.group_rows_by(mt[gene_field]).aggregate(
        hets=hl.agg.count_where(mt.GT.is_het()),
        homs=hl.agg.count_where(mt.GT.is_hom_var()),
        multi_het=hl.agg.count_where(mt.GT.is_het()) >= 2,
    )

    n_genes = mt_genes.count_rows()
    logger.info("  %d unique genes after aggregation", n_genes)

    return mt_genes


def compute_geneset_burden_mt(
    mt: hl.MatrixTable,
    gene_sets: Dict[str, List[str]],
    gene_field: str = "SYMBOL",
    genotype_aggregation: str = "hets",
    variant_filter: Optional["VariantFilter"] = None,
    max_af: Optional[float] = None,
    min_score: Optional[float] = None,
    consequences: Optional[List[str]] = None,
    normalize_by_length: bool = False,
    gene_lengths: Optional[Dict[str, float]] = None,
    min_gene_set_size: int = 0,
) -> Optional["hl.MatrixTable"]:
    """Compute per-sample burden for each gene set, returning a MatrixTable.

    This follows the two-step aggregation pattern:
    1. First aggregate variants → genes (per sample)
    2. Then aggregate genes → gene sets (per sample)

    The input MatrixTable is expected to be pre-filtered (QC, sample/variant
    filtering done upstream).  An optional ``variant_filter`` can be supplied
    for variant-class selection (e.g., restrict to LoF consequences).

    Parameters
    ----------
    mt : hl.MatrixTable
        Annotated variant MatrixTable with genotypes.  Expected to be
        pre-filtered for QC; only variant-class selection is applied here.
    gene_sets : Dict[str, List[str]]
        Dictionary mapping gene set names to gene lists.
    gene_field : str
        Row field containing gene symbol.
    genotype_aggregation : str
        One of "hets", "homs", "multi_het", "homs_multi_het".
        Deprecated aliases "chets" and "homs_chets" are accepted with a warning.
        "multi_het" counts genes with >= 2 heterozygous qualifying variants
        (approximates compound-het; over-counts when variants are in cis).
    variant_filter : VariantFilter, optional
        Variant filter criteria for variant-class selection.  When provided,
        legacy parameters (max_af, min_score, consequences) are ignored.
    max_af : float, optional
        Deprecated. Use ``variant_filter`` instead.
    min_score : float, optional
        Deprecated. Use ``variant_filter`` instead.
    consequences : List[str], optional
        Deprecated. Use ``variant_filter`` instead.
    normalize_by_length : bool
        If True, normalize per-gene burden by gene length before
        aggregating into gene sets.  This controls for the tendency of
        longer genes to accumulate more rare variants.  Burden becomes
        a float (rate per kb) rather than integer (count).
    gene_lengths : Dict[str, float], optional
        Gene symbol to CDS length in bases.  When *None* and
        ``normalize_by_length=True``, uses qualifying variant site count
        per gene as proxy.
    min_gene_set_size : int
        Minimum number of genes required per gene set.  Gene sets smaller
        than this are dropped before burden computation.  Default 0 (no
        filtering).

    Returns
    -------
    Optional[hl.MatrixTable]
        MatrixTable with rows = gene sets, cols = samples, entry = burden.
        Row fields include ``gene_set_size`` (original size),
        ``n_genes_found`` (genes with qualifying variants in the cohort),
        and ``gene_coverage_pct`` (``n_genes_found / gene_set_size * 100``).
        Returns ``None`` when no qualifying variants are found in gene set
        genes.

    Examples
    --------
    >>> gene_sets = {
    ...     "microglia": ["TREM2", "CD33", "ABI3"],
    ...     "astrocytes": ["GFAP", "AQP4", "S100B"]
    ... }
    >>> # No filtering (MT is pre-filtered):
    >>> mt_burden = compute_geneset_burden_mt(mt, gene_sets)
    >>> # With variant-class selection:
    >>> vf = VariantFilter(consequences=["stop_gained", "frameshift"])
    >>> mt_burden = compute_geneset_burden_mt(mt, gene_sets, variant_filter=vf)
    """
    _require_hail()

    # Resolve deprecated aliases
    genotype_aggregation = _resolve_genotype_aggregation(genotype_aggregation)

    # Pre-filter gene sets by minimum size
    if min_gene_set_size > 0:
        n_before = len(gene_sets)
        gene_sets = {k: v for k, v in gene_sets.items() if len(v) >= min_gene_set_size}
        n_dropped = n_before - len(gene_sets)
        if n_dropped > 0:
            logger.info(
                "Filtered %d/%d gene sets with fewer than %d genes",
                n_dropped,
                n_before,
                min_gene_set_size,
            )
        if not gene_sets:
            logger.warning(
                "All gene sets filtered out by min_gene_set_size=%d",
                min_gene_set_size,
            )
            return None

    gs_sizes = sorted(len(g) for g in gene_sets.values())
    n_gs = len(gene_sets)
    logger.info("Computing burden for %d gene sets", n_gs)
    if n_gs > 0:
        median_size = gs_sizes[n_gs // 2]
        logger.info(
            "  Gene set sizes: min=%d, median=%d, max=%d",
            gs_sizes[0],
            median_size,
            gs_sizes[-1],
        )
    logger.info(f"  Genotype aggregation: {genotype_aggregation}")

    # Validate genotype aggregation method
    valid_methods = ["hets", "homs", "multi_het", "homs_multi_het"]
    if genotype_aggregation not in valid_methods:
        raise ValueError(
            f"Invalid genotype_aggregation: {genotype_aggregation}. "
            f"Must be one of: {valid_methods}"
        )

    # Check required field exists
    if gene_field not in mt.row:
        raise ValueError(f"Gene field '{gene_field}' not found in MatrixTable")

    # Resolve variant filtering
    _has_legacy_params = (
        max_af is not None or min_score is not None or consequences is not None
    )

    if variant_filter is not None:
        if _has_legacy_params:
            logger.warning(
                "Both variant_filter and legacy parameters (max_af/min_score/consequences) "
                "provided; variant_filter takes precedence."
            )
        logger.info("Applying variant filter:")
        filter_expr = variant_filter.to_hail_expr(mt)
        mt = mt.filter_rows(filter_expr)
        n_variants = mt.count_rows()
        logger.info(f"  {n_variants} qualifying variants after filtering")
    elif _has_legacy_params:
        warnings.warn(
            "The max_af, min_score, and consequences parameters are deprecated. "
            "Use variant_filter=VariantFilter(...) instead, or pre-filter the "
            "MatrixTable before calling compute_geneset_burden_mt().",
            DeprecationWarning,
            stacklevel=2,
        )
        legacy_vf = VariantFilter(
            max_af=max_af if max_af is not None else 1.0,
            min_score=min_score,
            consequences=consequences,
            pass_only=False,
            min_gq=0,
            min_dp=0,
        )
        filter_expr = legacy_vf.to_hail_expr(mt)
        mt = mt.filter_rows(filter_expr)
        n_variants = mt.count_rows()
        logger.info(f"  {n_variants} qualifying variants after legacy filtering")
    else:
        logger.info("No variant filtering applied (MT assumed pre-filtered)")

    # Gene-length normalization setup
    _gene_length_ht = None
    if normalize_by_length:
        if gene_lengths is not None:
            logger.info(
                "Gene-length normalization: using provided CDS lengths (%d genes)",
                len(gene_lengths),
            )
            _gene_length_ht = hl.Table.parallelize(
                [
                    hl.struct(gene=k, _gene_length=float(v) / 1000.0)
                    for k, v in gene_lengths.items()
                ],
                schema=hl.tstruct(gene=hl.tstr, _gene_length=hl.tfloat64),
            ).key_by("gene")
        else:
            logger.info(
                "Gene-length normalization: using qualifying variant site "
                "count per gene as proxy (provide gene_lengths for "
                "CDS-based normalization)"
            )
            _rows = mt.rows()
            _gene_length_ht = _rows.group_by(gene=_rows[gene_field]).aggregate(
                _gene_length=hl.float64(hl.agg.count())
            )

    # Create gene → gene_sets mapping
    logger.info("Creating gene to gene set mapping...")
    gene_to_sets = {}
    for gs_name, genes in gene_sets.items():
        for gene in genes:
            if gene not in gene_to_sets:
                gene_to_sets[gene] = []
            gene_to_sets[gene].append(gs_name)

    logger.info(f"  {len(gene_to_sets)} unique genes across all gene sets")

    # Convert to Hail Table
    gene_to_sets_ht = hl.Table.parallelize(
        [hl.struct(gene=k, gene_set_ids=v) for k, v in gene_to_sets.items()],
        schema=hl.tstruct(gene=hl.tstr, gene_set_ids=hl.tarray(hl.tstr)),
    ).key_by("gene")

    # Annotate MT with gene set membership
    mt = mt.annotate_rows(gene_set_ids=gene_to_sets_ht[mt[gene_field]].gene_set_ids)

    # Filter to only genes in our gene sets
    mt = mt.filter_rows(hl.is_defined(mt.gene_set_ids))
    n_variants_in_sets = mt.count_rows()
    logger.info(f"  {n_variants_in_sets} variants in gene set genes")

    if n_variants_in_sets == 0:
        logger.warning(
            "No variants found in gene set genes after filtering. "
            "Check that gene field '%s' and gene set IDs match.",
            gene_field,
        )
        return None

    # STEP 1: Aggregate variants → genes per sample
    logger.info("Step 1: Aggregating variants to genes per sample...")
    mt_genes = mt.group_rows_by(mt[gene_field]).aggregate(
        hets=hl.agg.count_where(mt.GT.is_het()),
        homs=hl.agg.count_where(mt.GT.is_hom_var()),
        # multi_het: genes with >= 2 heterozygous qualifying variants.
        # This approximates compound-het but over-counts when variants are in cis.
        multi_het=hl.agg.count_where(mt.GT.is_het()) >= 2,
    )

    # Annotate with gene length (for normalization)
    if normalize_by_length and _gene_length_ht is not None:
        mt_genes = mt_genes.annotate_rows(
            _gene_length=_gene_length_ht[mt_genes[gene_field]]._gene_length
        )
        # Default to 1.0 for genes without length data
        mt_genes = mt_genes.annotate_rows(
            _gene_length=hl.if_else(
                hl.is_defined(mt_genes._gene_length) & (mt_genes._gene_length > 0),
                mt_genes._gene_length,
                1.0,
            )
        )

    # Re-annotate with gene set membership (lost in grouping)
    mt_genes = mt_genes.annotate_rows(
        gene_set_ids=gene_to_sets_ht[mt_genes[gene_field]].gene_set_ids
    )

    # Explode by gene set (so each gene can contribute to multiple sets)
    mt_genes = mt_genes.annotate_rows(
        gene_set_ids=hl.or_else(mt_genes.gene_set_ids, hl.empty_array(hl.tstr))
    )
    mt_genes = mt_genes.explode_rows(mt_genes.gene_set_ids)

    # STEP 2: Aggregate genes → gene sets per sample
    logger.info("Step 2: Aggregating genes to gene sets per sample...")

    # Select aggregation method
    if normalize_by_length:
        # Rate-based aggregation: divide by gene length
        if genotype_aggregation == "hets":
            agg_expr = hl.agg.sum(
                hl.if_else(mt_genes.hets > 0, 1.0 / mt_genes._gene_length, 0.0)
            )
        elif genotype_aggregation == "homs":
            agg_expr = hl.agg.sum(
                hl.if_else(mt_genes.homs > 0, 1.0 / mt_genes._gene_length, 0.0)
            )
        elif genotype_aggregation == "multi_het":
            agg_expr = hl.agg.sum(
                hl.if_else(mt_genes.multi_het, 1.0 / mt_genes._gene_length, 0.0)
            )
        elif genotype_aggregation == "homs_multi_het":
            agg_expr = hl.agg.sum(
                hl.if_else(
                    mt_genes.multi_het | (mt_genes.homs > 0),
                    1.0 / mt_genes._gene_length,
                    0.0,
                )
            )
    elif genotype_aggregation == "hets":
        agg_expr = hl.int(hl.agg.sum(hl.if_else(mt_genes.hets > 0, 1, 0)))
    elif genotype_aggregation == "homs":
        agg_expr = hl.int(hl.agg.sum(hl.if_else(mt_genes.homs > 0, 1, 0)))
    elif genotype_aggregation == "multi_het":
        agg_expr = hl.int(hl.agg.sum(hl.if_else(mt_genes.multi_het, 1, 0)))
    elif genotype_aggregation == "homs_multi_het":
        agg_expr = hl.int(
            hl.agg.sum(hl.if_else(mt_genes.multi_het | (mt_genes.homs > 0), 1, 0))
        )

    # Count genes found per gene set (before grouping collapses gene info).
    # Must use the rows Table's own field reference — not mt_genes — to
    # avoid "Cannot combine expressions from different source objects".
    _rows_ht = mt_genes.rows()
    _genes_per_set_ht = _rows_ht.group_by(
        gene_set_name=_rows_ht.gene_set_ids
    ).aggregate(n_genes_found=hl.agg.count())

    mt_burden = mt_genes.group_rows_by(gene_set_name=mt_genes.gene_set_ids).aggregate(
        burden=agg_expr
    )

    # Annotate with gene set metrics
    _gene_set_size_ht = hl.Table.parallelize(
        [
            hl.struct(gene_set_name=name, gene_set_size=len(genes))
            for name, genes in gene_sets.items()
        ],
        schema=hl.tstruct(gene_set_name=hl.tstr, gene_set_size=hl.tint32),
    ).key_by("gene_set_name")

    mt_burden = mt_burden.annotate_rows(
        gene_set_size=_gene_set_size_ht[mt_burden.gene_set_name].gene_set_size,
        n_genes_found=hl.int32(
            _genes_per_set_ht[mt_burden.gene_set_name].n_genes_found
        ),
    )
    mt_burden = mt_burden.annotate_rows(
        gene_coverage_pct=hl.format(
            "%.1f",
            hl.float64(mt_burden.n_genes_found) / mt_burden.gene_set_size * 100,
        )
    )

    n_gene_sets_final = mt_burden.count_rows()
    n_samples = mt_burden.count_cols()
    logger.info(
        f"Burden matrix created: {n_gene_sets_final} gene sets × {n_samples} samples"
    )

    # Warn about low-coverage gene sets
    n_low_cov = mt_burden.filter_rows(
        mt_burden.n_genes_found < (mt_burden.gene_set_size / 2)
    ).count_rows()
    if n_low_cov > 0:
        logger.warning(
            "%d/%d gene sets have <50%% gene coverage in the cohort — "
            "results may be unreliable for these sets",
            n_low_cov,
            n_gene_sets_final,
        )

    # Check for gene sets with zero burden across all samples
    _burden_sums = mt_burden.annotate_rows(_total_burden=hl.agg.sum(mt_burden.burden))
    n_zero = _burden_sums.filter_rows(_burden_sums._total_burden == 0).count_rows()
    if n_zero > 0:
        logger.warning(
            "%d/%d gene sets have zero burden across all samples",
            n_zero,
            n_gene_sets_final,
        )

    return mt_burden


def logistic_burden_test(
    mt_burden: hl.MatrixTable,
    phenotype_field: str,
    covariates: Optional[List[str]] = None,
    pass_through: Optional[List[str]] = None,
) -> hl.Table:
    """Run Hail-native logistic regression for burden testing.

    Uses hl.logistic_regression_rows() internally, following the pattern
    from alphagenome/logreg_burden_test.py and hvantk/utils/stats.py.

    Parameters
    ----------
    mt_burden : hl.MatrixTable
        Burden MatrixTable with rows = gene sets, cols = samples, entry = burden
    phenotype_field : str
        Column field containing binary phenotype (True/False or 1/0)
    covariates : List[str], optional
        Column fields to use as covariates (e.g., ["PC1", "PC2", "sex"])
    pass_through : List[str], optional
        Row fields to include in output

    Returns
    -------
    hl.Table
        Results table with columns:
        - gene_set_name
        - beta (coefficient)
        - standard_error
        - z_stat
        - p_value
        - odds_ratio (exp(beta))
        - ci_lower, ci_upper (95% CI)
        - fit (convergence info)

    Examples
    --------
    >>> result_ht = logistic_burden_test(
    ...     mt_burden,
    ...     phenotype_field="is_case",
    ...     covariates=["PC1", "PC2", "PC3", "sex"]
    ... )
    >>> result_ht.show()
    """
    _require_hail()
    n_gene_sets = mt_burden.count_rows()
    logger.info("Running Hail-native logistic regression (%d gene sets)", n_gene_sets)
    logger.info(f"  Phenotype: {phenotype_field}")

    # Check phenotype field exists
    if phenotype_field not in mt_burden.col:
        raise ValueError(f"Phenotype field '{phenotype_field}' not found in columns")

    # Build covariate expressions
    cov_exprs = [1.0]  # Intercept
    if covariates:
        logger.info(f"  Covariates: {', '.join(covariates)}")
        for cov in covariates:
            if cov not in mt_burden.col:
                raise ValueError(f"Covariate '{cov}' not found in columns")
            cov_exprs.append(mt_burden[cov])
    else:
        logger.info("  No covariates")

    # Convert phenotype to float if needed
    pheno_expr = mt_burden[phenotype_field]
    if pheno_expr.dtype != hl.tfloat64 and pheno_expr.dtype != hl.tfloat32:
        pheno_expr = hl.float(pheno_expr)

    # Run logistic regression
    logger.info("Running regression...")
    result = hl.logistic_regression_rows(
        test="wald",
        y=pheno_expr,
        x=hl.float(mt_burden.burden),
        covariates=cov_exprs,
        pass_through=pass_through or [],
    )

    # Annotate with additional statistics
    result = result.annotate(
        odds_ratio=hl.exp(result.beta),
        ci_lower=hl.exp(result.beta - 1.96 * result.standard_error),
        ci_upper=hl.exp(result.beta + 1.96 * result.standard_error),
    )

    n_tested = result.count()
    logger.info("Regression complete: %d gene sets tested", n_tested)

    # Check for convergence issues (NaN p-values)
    n_nan_p = result.filter(hl.is_nan(result.p_value)).count()
    if n_nan_p > 0:
        logger.warning(
            "Logistic regression produced NaN p-values for %d/%d gene sets "
            "(possible convergence failure)",
            n_nan_p,
            n_tested,
        )

    # Report nominally significant results
    n_nominal = result.filter(result.p_value < 0.05).count()
    logger.info("  %d gene sets nominally significant (p < 0.05)", n_nominal)

    return result


def linear_burden_test(
    mt_burden: hl.MatrixTable,
    phenotype_field: str,
    covariates: Optional[List[str]] = None,
    pass_through: Optional[List[str]] = None,
) -> hl.Table:
    """Run Hail-native linear regression for burden testing.

    Uses hl.linear_regression_rows() for continuous phenotypes.

    Parameters
    ----------
    mt_burden : hl.MatrixTable
        Burden MatrixTable with rows = gene sets, cols = samples, entry = burden
    phenotype_field : str
        Column field containing continuous phenotype
    covariates : List[str], optional
        Column fields to use as covariates
    pass_through : List[str], optional
        Row fields to include in output

    Returns
    -------
    hl.Table
        Results table with columns:
        - gene_set_name
        - beta (coefficient)
        - standard_error
        - t_stat
        - p_value

    Examples
    --------
    >>> result_ht = linear_burden_test(
    ...     mt_burden,
    ...     phenotype_field="cognitive_score",
    ...     covariates=["age", "sex", "PC1", "PC2"]
    ... )
    """
    _require_hail()
    n_gene_sets = mt_burden.count_rows()
    logger.info("Running Hail-native linear regression (%d gene sets)", n_gene_sets)
    logger.info(f"  Phenotype: {phenotype_field}")

    # Check phenotype field exists
    if phenotype_field not in mt_burden.col:
        raise ValueError(f"Phenotype field '{phenotype_field}' not found in columns")

    # Build covariate expressions
    cov_exprs = [1.0]  # Intercept
    if covariates:
        logger.info(f"  Covariates: {', '.join(covariates)}")
        for cov in covariates:
            if cov not in mt_burden.col:
                raise ValueError(f"Covariate '{cov}' not found in columns")
            cov_exprs.append(mt_burden[cov])
    else:
        logger.info("  No covariates")

    # Convert phenotype to float if needed
    pheno_expr = mt_burden[phenotype_field]
    if pheno_expr.dtype != hl.tfloat64 and pheno_expr.dtype != hl.tfloat32:
        pheno_expr = hl.float(pheno_expr)

    # Run linear regression
    logger.info("Running regression...")
    result = hl.linear_regression_rows(
        y=pheno_expr,
        x=hl.float(mt_burden.burden),
        covariates=cov_exprs,
        pass_through=pass_through or [],
    )

    n_tested = result.count()
    logger.info("Regression complete: %d gene sets tested", n_tested)

    # Check for convergence issues (NaN p-values)
    n_nan_p = result.filter(hl.is_nan(result.p_value)).count()
    if n_nan_p > 0:
        logger.warning(
            "Linear regression produced NaN p-values for %d/%d gene sets "
            "(possible convergence failure)",
            n_nan_p,
            n_tested,
        )

    # Report nominally significant results
    n_nominal = result.filter(result.p_value < 0.05).count()
    logger.info("  %d gene sets nominally significant (p < 0.05)", n_nominal)

    return result


@algorithm(name="burden_analysis", backends=[Backend.HAIL])
def run_burden_analysis(
    cohort_mt: hl.MatrixTable,
    gene_sets: Dict[str, List[str]],
    phenotype_ht: hl.Table,
    phenotype_field: str = "is_case",
    covariate_fields: Optional[List[str]] = None,
    phenotype_type: str = "binary",
    gene_field: str = "SYMBOL",
    genotype_aggregation: str = "hets",
    variant_filter: Optional["VariantFilter"] = None,
    max_af: Optional[float] = None,
    min_score: Optional[float] = None,
    consequences: Optional[List[str]] = None,
    normalize_by_length: bool = False,
    gene_lengths: Optional[Dict[str, float]] = None,
    min_carriers: int = 0,
    min_gene_set_size: int = 0,
) -> Optional["hl.Table"]:
    """Complete burden analysis pipeline.

    This is the main entry point that orchestrates:
    1. Burden computation per gene set
    2. Phenotype/covariate annotation
    3. Regression testing

    The cohort MatrixTable is expected to be pre-filtered (QC, sample/variant
    filtering done upstream).  An optional ``variant_filter`` can be supplied
    for variant-class selection.

    Parameters
    ----------
    cohort_mt : hl.MatrixTable
        Pre-filtered cohort MatrixTable with genotypes and gene annotation.
    gene_sets : Dict[str, List[str]]
        Gene sets to test.
    phenotype_ht : hl.Table
        Table with sample phenotypes, keyed by sample ID.
    phenotype_field : str
        Field containing phenotype.
    covariate_fields : List[str], optional
        Fields to use as covariates.
    phenotype_type : str
        "binary" or "continuous".
    gene_field : str
        Row field for gene symbol in cohort_mt.
    genotype_aggregation : str
        Genotype aggregation method.
    variant_filter : VariantFilter, optional
        Variant filter for variant-class selection.
    max_af : float, optional
        Deprecated. Use ``variant_filter`` instead.
    min_score : float, optional
        Deprecated. Use ``variant_filter`` instead.
    consequences : List[str], optional
        Deprecated. Use ``variant_filter`` instead.
    normalize_by_length : bool
        Normalize per-gene burden by gene length before aggregation.
    gene_lengths : Dict[str, float], optional
        Gene symbol to CDS length in bases for normalization.
    min_carriers : int
        Minimum number of samples carrying qualifying variants for a gene
        set to be included in regression.  Gene sets with fewer carriers
        are filtered out (regression would be uninformative).  Default 0
        (no filtering).
    min_gene_set_size : int
        Minimum number of genes required per gene set.  Gene sets smaller
        than this are dropped before burden computation.  Default 0.

    Returns
    -------
    Optional[hl.Table]
        Results with p-values, odds ratios, etc. for each gene set.
        Includes ``gene_set_size``, ``n_genes_found``, and
        ``gene_coverage_pct`` columns.
        Returns ``None`` when no testable gene sets remain (e.g., no
        qualifying variants, no phenotype overlap, or all gene sets
        filtered by ``min_carriers``).

    Examples
    --------
    >>> gene_sets = {
    ...     "microglia": ["TREM2", "CD33", "ABI3"],
    ...     "astrocytes": ["GFAP", "AQP4", "S100B"]
    ... }
    >>> pheno_ht = hl.import_table("phenotypes.tsv", impute=True).key_by("sample_id")
    >>> results = run_burden_analysis(
    ...     cohort_mt=mt,
    ...     gene_sets=gene_sets,
    ...     phenotype_ht=pheno_ht,
    ...     phenotype_field="is_case",
    ...     covariate_fields=["PC1", "PC2", "sex"]
    ... )
    """
    _require_hail()
    logger.info("=" * 60)
    logger.info("BURDEN ANALYSIS PIPELINE")
    logger.info("=" * 60)

    # Validate inputs
    if phenotype_type not in ["binary", "continuous"]:
        raise ValueError(f"Invalid phenotype_type: {phenotype_type}")

    # Step 1: Compute burden matrix
    logger.info("\n[1/3] Computing burden matrix...")
    t_burden = time.time()
    mt_burden = compute_geneset_burden_mt(
        mt=cohort_mt,
        gene_sets=gene_sets,
        gene_field=gene_field,
        genotype_aggregation=genotype_aggregation,
        variant_filter=variant_filter,
        max_af=max_af,
        min_score=min_score,
        consequences=consequences,
        normalize_by_length=normalize_by_length,
        gene_lengths=gene_lengths,
        min_gene_set_size=min_gene_set_size,
    )
    logger.info("  Burden computation took %.1fs", time.time() - t_burden)

    if mt_burden is None:
        logger.warning("Burden computation returned no results — skipping regression.")
        return None

    # Step 2: Annotate with phenotypes and covariates
    logger.info("\n[2/3] Annotating with phenotypes and covariates...")
    logger.info(f"  Joining phenotype table (keyed by: {phenotype_ht.key})")

    # Get the column key from the burden MT (may not be 's' after grouping)
    col_key = list(mt_burden.col_key)[0] if mt_burden.col_key else None
    if col_key is None:
        raise ValueError("Burden MatrixTable has no column key")

    mt_burden = mt_burden.annotate_cols(**phenotype_ht[mt_burden[col_key]])

    # Filter to samples with phenotype
    mt_burden = mt_burden.filter_cols(hl.is_defined(mt_burden[phenotype_field]))

    n_samples_with_pheno = mt_burden.count_cols()
    logger.info(f"  {n_samples_with_pheno} samples with phenotype data")

    if n_samples_with_pheno == 0:
        logger.warning(
            "No samples with phenotype data after joining phenotype table. "
            "Check that sample IDs in the phenotype table match the cohort."
        )
        return None

    # Filter gene sets by min_carriers
    if min_carriers > 0:
        mt_burden = mt_burden.annotate_rows(
            _n_carriers=hl.agg.count_where(mt_burden.burden > 0)
        )
        n_before = mt_burden.count_rows()
        mt_burden = mt_burden.filter_rows(mt_burden._n_carriers >= min_carriers)
        mt_burden = mt_burden.drop("_n_carriers")
        n_after = mt_burden.count_rows()
        n_dropped = n_before - n_after
        if n_dropped > 0:
            logger.info(
                "  Filtered %d/%d gene sets with fewer than %d carriers",
                n_dropped,
                n_before,
                min_carriers,
            )
        if n_after == 0:
            logger.warning(
                "All gene sets filtered out by min_carriers=%d — no testable gene sets remain.",
                min_carriers,
            )
            return None

    # Step 3: Run regression
    logger.info(f"\n[3/3] Running {phenotype_type} regression...")
    t_regression = time.time()

    _pass_through = ["gene_set_size", "n_genes_found", "gene_coverage_pct"]

    if phenotype_type == "binary":
        result = logistic_burden_test(
            mt_burden,
            phenotype_field=phenotype_field,
            covariates=covariate_fields,
            pass_through=_pass_through,
        )
    else:
        result = linear_burden_test(
            mt_burden,
            phenotype_field=phenotype_field,
            covariates=covariate_fields,
            pass_through=_pass_through,
        )
    logger.info("  Regression took %.1fs", time.time() - t_regression)

    logger.info("\n" + "=" * 60)
    logger.info("BURDEN ANALYSIS COMPLETE")
    logger.info("=" * 60)

    return result


def build_variant_classes_from_presets(
    class_names: List[str],
    base_filter: Optional["VariantFilter"] = None,
) -> Dict[str, "VariantFilter"]:
    """Build variant class filters from preset or custom names.

    Each preset defines the consequence filter (and optionally a score
    threshold) for a variant class.  Non-consequence settings (max_af,
    field names, QC thresholds) are inherited from *base_filter* when
    provided, or default to permissive values (no filtering) since the
    MatrixTable is expected to be pre-filtered.

    Names not found in ``VARIANT_CLASS_PRESETS`` are treated as custom
    classes: the name itself is used as the consequence value.  This
    allows filtering on pre-computed consequence group fields (e.g.,
    ``csq_group`` with values ``hcLOF``, ``missC``, ``syn``).

    Parameters
    ----------
    class_names : List[str]
        Names of variant classes.  Recognised preset names:
        ``"lof"``, ``"missense_constrained"``, ``"synonymous"``.
        Any other name is treated as a custom class whose single
        qualifying consequence equals the name itself.
    base_filter : VariantFilter, optional
        Base filter whose non-consequence settings are inherited.  If
        *None*, uses permissive defaults (``max_af=1.0``,
        ``pass_only=False``, ``min_gq=0``, ``min_dp=0``).

    Returns
    -------
    Dict[str, VariantFilter]
        Mapping of class name to configured VariantFilter.

    Examples
    --------
    >>> # Using presets (VEP consequence terms)
    >>> classes = build_variant_classes_from_presets(
    ...     ["lof", "missense_constrained", "synonymous"]
    ... )
    >>> classes["lof"].consequences
    ['stop_gained', 'frameshift_variant', 'splice_donor_variant', 'splice_acceptor_variant']

    >>> # Using custom names (pre-computed csq_group values)
    >>> from hvantk.algorithms.enrichex.burden import VariantFilter
    >>> base = VariantFilter(consequence_field="csq_group", af_field="internal_af")
    >>> classes = build_variant_classes_from_presets(
    ...     ["hcLOF", "missC", "syn"], base_filter=base
    ... )
    >>> classes["hcLOF"].consequences
    ['hcLOF']
    """
    result: Dict[str, VariantFilter] = {}
    for name in class_names:
        preset = VARIANT_CLASS_PRESETS.get(name)

        if preset is not None:
            # Known preset — use its consequence list and optional score
            consequences = preset.get("consequences")
            min_score = preset.get(
                "min_score", base_filter.min_score if base_filter else None
            )
        else:
            # Custom class — name IS the consequence value
            logger.info(
                "Variant class '%s' is not a preset; treating as custom "
                "class with consequence=['%s']",
                name,
                name,
            )
            consequences = [name]
            min_score = base_filter.min_score if base_filter else None

        if base_filter is not None:
            vf = VariantFilter(
                max_af=base_filter.max_af,
                min_score=min_score,
                consequences=consequences,
                pass_only=base_filter.pass_only,
                min_gq=base_filter.min_gq,
                min_dp=base_filter.min_dp,
                af_field=base_filter.af_field,
                score_field=base_filter.score_field,
                consequence_field=base_filter.consequence_field,
            )
        else:
            vf = VariantFilter(
                max_af=1.0,
                min_score=min_score,
                consequences=consequences,
                pass_only=False,
                min_gq=0,
                min_dp=0,
            )

        result[name] = vf

    return result


def _check_synonymous_control(
    synonymous_result: "hl.Table", alpha: float = 0.05
) -> None:
    """Check synonymous results for potential confounding.

    Logs a warning if any gene set shows nominal significance in the
    synonymous variant class, which may indicate gene-length bias or
    other systematic confounding.

    Parameters
    ----------
    synonymous_result : hl.Table
        Burden regression results for the synonymous variant class.
    alpha : float
        Nominal significance threshold (default 0.05).
    """
    _require_hail()
    sig = synonymous_result.filter(synonymous_result.p_value < alpha)
    n_sig = sig.count()
    if n_sig > 0:
        gene_sets = sig.gene_set_name.collect()
        logger.warning(
            "SYNONYMOUS CONTROL CHECK: %d gene set(s) show nominal "
            "significance (p < %s) for synonymous variants: %s. "
            "This may indicate gene-length bias or other confounding.",
            n_sig,
            alpha,
            gene_sets,
        )


@algorithm(name="stratified_burden_analysis", backends=[Backend.HAIL])
def run_stratified_burden_analysis(
    cohort_mt: "hl.MatrixTable",
    gene_sets: Dict[str, List[str]],
    phenotype_ht: "hl.Table",
    variant_classes: Dict[str, "VariantFilter"],
    phenotype_field: str = "is_case",
    covariate_fields: Optional[List[str]] = None,
    phenotype_type: str = "binary",
    gene_field: str = "SYMBOL",
    genotype_aggregation: str = "hets",
    normalize_by_length: bool = False,
    gene_lengths: Optional[Dict[str, float]] = None,
    min_carriers: int = 0,
    min_gene_set_size: int = 0,
) -> Dict[str, "hl.Table"]:
    """Run burden analysis stratified by variant class.

    For each variant class (e.g., LoF, constrained missense, synonymous),
    runs an independent burden analysis against the same gene sets using
    the class-specific ``VariantFilter``.  The synonymous class serves as
    a negative control: a warning is logged if any gene set shows nominal
    significance for synonymous variants.

    Parameters
    ----------
    cohort_mt : hl.MatrixTable
        Pre-filtered cohort MatrixTable with genotypes and gene annotation.
    gene_sets : Dict[str, List[str]]
        Gene sets to test.
    phenotype_ht : hl.Table
        Table with sample phenotypes, keyed by sample ID.
    variant_classes : Dict[str, VariantFilter]
        Mapping of variant class name to filter criteria.  Use
        ``build_variant_classes_from_presets()`` for convenience or
        construct manually.
    phenotype_field : str
        Field containing phenotype.
    covariate_fields : List[str], optional
        Fields to use as covariates.
    phenotype_type : str
        ``"binary"`` or ``"continuous"``.
    gene_field : str
        Row field for gene symbol in cohort_mt.
    genotype_aggregation : str
        Genotype aggregation method.
    normalize_by_length : bool
        Normalize per-gene burden by gene length before aggregation.
    gene_lengths : Dict[str, float], optional
        Gene symbol to CDS length in bases for normalization.
    min_carriers : int
        Minimum number of carriers per gene set for inclusion in
        regression.  Passed through to ``run_burden_analysis()``.

    Returns
    -------
    Dict[str, hl.Table]
        Mapping of variant class name to results table.  Each table has
        an additional ``variant_class`` field.  Classes where no qualifying
        variants were found or no testable gene sets remain are omitted
        with a logged warning.

    Examples
    --------
    >>> from hvantk.algorithms.enrichex.burden import (
    ...     VariantFilter,
    ...     run_stratified_burden_analysis,
    ... )
    >>> variant_classes = {
    ...     "lof": VariantFilter(
    ...         consequences=["stop_gained", "frameshift_variant"],
    ...         pass_only=False, min_gq=0, min_dp=0,
    ...     ),
    ...     "synonymous": VariantFilter(
    ...         consequences=["synonymous_variant"],
    ...         pass_only=False, min_gq=0, min_dp=0,
    ...     ),
    ... }
    >>> results = run_stratified_burden_analysis(
    ...     cohort_mt=mt,
    ...     gene_sets=gene_sets,
    ...     phenotype_ht=pheno_ht,
    ...     variant_classes=variant_classes,
    ... )
    >>> for cls, ht in results.items():
    ...     print(cls, ht.count())
    """
    _require_hail()

    logger.info("=" * 60)
    logger.info("STRATIFIED BURDEN ANALYSIS")
    logger.info("=" * 60)
    logger.info("Variant classes: %s", list(variant_classes.keys()))

    results: Dict[str, Any] = {}

    for class_name, vf in variant_classes.items():
        logger.info("\n" + "-" * 60)
        logger.info("VARIANT CLASS: %s", class_name)
        logger.info("-" * 60)

        result = run_burden_analysis(
            cohort_mt=cohort_mt,
            gene_sets=gene_sets,
            phenotype_ht=phenotype_ht,
            variant_filter=vf,
            phenotype_field=phenotype_field,
            covariate_fields=covariate_fields,
            phenotype_type=phenotype_type,
            gene_field=gene_field,
            genotype_aggregation=genotype_aggregation,
            normalize_by_length=normalize_by_length,
            gene_lengths=gene_lengths,
            min_carriers=min_carriers,
            min_gene_set_size=min_gene_set_size,
        )
        if result is not None:
            result = result.annotate(variant_class=class_name)
            results[class_name] = result
        else:
            logger.warning(
                "Skipping variant class '%s': no testable results.", class_name
            )

    # Synonymous negative control check
    if "synonymous" in results:
        _check_synonymous_control(results["synonymous"])

    logger.info("\n" + "=" * 60)
    logger.info(
        "STRATIFIED BURDEN ANALYSIS COMPLETE — %d/%d classes with results",
        len(results),
        len(variant_classes),
    )
    logger.info("=" * 60)

    return results


# ---------------------------------------------------------------------------
# Competitive gene-set testing (permutation-based)
# ---------------------------------------------------------------------------


def _logistic_z_statistic(X: "np.ndarray", y: "np.ndarray") -> float:
    """Compute z-statistic for the first predictor via logistic IRLS.

    Parameters
    ----------
    X : np.ndarray, shape (n, p)
        Design matrix.  The first column is the predictor of interest
        (burden); remaining columns are intercept and covariates.
    y : np.ndarray, shape (n,)
        Binary outcome (0/1).

    Returns
    -------
    float
        z-statistic, or ``nan`` if regression fails to converge.
    """
    import numpy as np
    from scipy.special import expit

    _n, p = X.shape
    beta = np.zeros(p)

    for _ in range(25):
        eta = X @ beta
        # Clip to prevent overflow in expit
        eta = np.clip(eta, -20, 20)
        mu = expit(eta)
        w = np.maximum(mu * (1 - mu), 1e-10)
        z = eta + (y - mu) / w
        Xw = X * np.sqrt(w)[:, None]
        XWX = Xw.T @ Xw
        try:
            beta_new = np.linalg.solve(XWX, (X * w[:, None]).T @ z)
        except np.linalg.LinAlgError:
            return float("nan")
        if np.max(np.abs(beta_new - beta)) < 1e-8:
            beta = beta_new
            break
        beta = beta_new

    mu = expit(np.clip(X @ beta, -20, 20))
    w = np.maximum(mu * (1 - mu), 1e-10)
    Xw = X * np.sqrt(w)[:, None]
    try:
        cov = np.linalg.inv(Xw.T @ Xw)
    except np.linalg.LinAlgError:
        return float("nan")

    se = float(np.sqrt(max(cov[0, 0], 0)))
    if se < 1e-10:
        return float("nan")
    return float(beta[0] / se)


def _linear_t_statistic(X: "np.ndarray", y: "np.ndarray") -> float:
    """Compute t-statistic for the first predictor via OLS.

    Parameters
    ----------
    X : np.ndarray, shape (n, p)
        Design matrix.  First column is the predictor of interest.
    y : np.ndarray, shape (n,)
        Continuous outcome.

    Returns
    -------
    float
        t-statistic, or ``nan`` if regression fails.
    """
    import numpy as np

    n, p = X.shape
    try:
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
    except np.linalg.LinAlgError:
        return float("nan")

    resid = y - X @ beta
    sigma2 = float(np.sum(resid**2) / max(n - p, 1))
    try:
        cov = sigma2 * np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        return float("nan")

    se = float(np.sqrt(max(cov[0, 0], 0)))
    if se < 1e-10:
        return float("nan")
    return float(beta[0] / se)


def _compute_test_statistic(
    X: "np.ndarray", y: "np.ndarray", phenotype_type: str
) -> float:
    """Dispatch to logistic or linear regression test statistic."""
    if phenotype_type == "binary":
        return _logistic_z_statistic(X, y)
    return _linear_t_statistic(X, y)


def _sample_length_matched(
    target_indices: "np.ndarray",
    gene_bins: "np.ndarray",
    bin_to_indices: Dict[int, List[int]],
    rng: "np.random.Generator",
) -> "np.ndarray":
    """Sample random genes matched by length bin.

    For each gene in *target_indices*, samples a random gene from the
    same length bin.

    Parameters
    ----------
    target_indices : np.ndarray
        Gene indices of the target gene set.
    gene_bins : np.ndarray
        Bin assignment for every gene in the background.
    bin_to_indices : dict
        Mapping from bin ID to list of gene indices in that bin.
    rng : np.random.Generator
        Random number generator.

    Returns
    -------
    np.ndarray
        Indices of the sampled genes.
    """
    import numpy as np

    sampled = []
    for idx in target_indices:
        b = int(gene_bins[idx])
        candidates = bin_to_indices.get(b, [int(idx)])
        sampled.append(rng.choice(candidates))
    return np.asarray(sampled, dtype=np.intp)


def permutation_burden_test(
    cohort_mt: "hl.MatrixTable",
    gene_sets: Dict[str, List[str]],
    phenotype_ht: "hl.Table",
    phenotype_field: str = "is_case",
    covariate_fields: Optional[List[str]] = None,
    phenotype_type: str = "binary",
    gene_field: str = "SYMBOL",
    genotype_aggregation: str = "hets",
    variant_filter: Optional["VariantFilter"] = None,
    n_permutations: int = 10000,
    length_matched: bool = True,
    seed: Optional[int] = None,
    normalize_by_length: bool = False,
    gene_lengths: Optional[Dict[str, float]] = None,
) -> "pd.DataFrame":
    """Competitive gene-set burden test via gene-label permutation.

    Tests whether each gene set is **more** burdened than random gene
    sets of the same size (competitive null), as opposed to the
    self-contained null tested by ``run_burden_analysis()``.  This
    controls for the baseline rate of rare variation and gene set size.

    Algorithm
    ---------
    1. Compute per-gene qualifying indicator per sample (in Hail).
    2. Collect per-sample qualifying gene sets to local memory.
    3. For each target gene set, compute burden and regression
       z/t-statistic (observed).
    4. Repeat *n_permutations* times: sample a random gene set of the
       same size, compute burden and test statistic.
    5. Empirical p-value =
       ``(# permutations with |stat| >= |observed| + 1) / (N + 1)``.

    When ``length_matched=True`` (default), random gene sets are sampled
    from length-matched bins using the number of qualifying samples per
    gene as a proxy for gene length / variant density.

    Parameters
    ----------
    cohort_mt : hl.MatrixTable
        Pre-filtered cohort MatrixTable.
    gene_sets : Dict[str, List[str]]
        Gene sets to test.
    phenotype_ht : hl.Table
        Phenotype table, keyed by sample ID.
    phenotype_field : str
        Phenotype column name.
    covariate_fields : List[str], optional
        Covariate column names.
    phenotype_type : str
        ``"binary"`` or ``"continuous"``.
    gene_field : str
        Row field for gene symbol.
    genotype_aggregation : str
        Genotype aggregation method.
    variant_filter : VariantFilter, optional
        Variant-class filter.
    n_permutations : int
        Number of permutations (default 10 000).
    length_matched : bool
        Use length-matched gene sampling (default True).
    seed : int, optional
        Random seed for reproducibility.
    normalize_by_length : bool
        Normalize per-gene qualifying contribution by gene length.
    gene_lengths : Dict[str, float], optional
        Gene symbol to CDS length in bases for normalization.

    Returns
    -------
    pd.DataFrame
        One row per gene set with columns: ``gene_set_name``,
        ``n_genes_tested``, ``observed_statistic``,
        ``empirical_p_value``, ``n_permutations``, ``n_exceeded``.
    """
    _require_hail()
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix, lil_matrix

    logger.info("=" * 60)
    logger.info("PERMUTATION BURDEN TEST")
    logger.info("=" * 60)
    logger.info("  Gene sets: %d", len(gene_sets))
    logger.info("  Permutations: %d", n_permutations)
    logger.info("  Length matched: %s", length_matched)

    # Resolve deprecated aliases
    genotype_aggregation = _resolve_genotype_aggregation(genotype_aggregation)
    valid_methods = ["hets", "homs", "multi_het", "homs_multi_het"]
    if genotype_aggregation not in valid_methods:
        raise ValueError(
            f"Invalid genotype_aggregation: {genotype_aggregation}. "
            f"Must be one of: {valid_methods}"
        )

    # Step 1: Compute per-gene burden MT (all genes, not just gene set genes)
    logger.info("Computing per-gene burden matrix...")
    mt_genes = compute_per_gene_burden_mt(cohort_mt, gene_field, variant_filter)

    # Determine qualification expression
    if genotype_aggregation == "hets":
        qual_expr = mt_genes.hets > 0
    elif genotype_aggregation == "homs":
        qual_expr = mt_genes.homs > 0
    elif genotype_aggregation == "multi_het":
        qual_expr = mt_genes.multi_het
    else:  # homs_multi_het
        qual_expr = mt_genes.multi_het | (mt_genes.homs > 0)

    mt_genes = mt_genes.annotate_entries(_qualifies=qual_expr)

    # Step 2: Annotate with phenotype and collect to local
    logger.info("Annotating with phenotype data...")
    col_key_name = list(mt_genes.col_key)[0]
    mt_genes = mt_genes.annotate_cols(**phenotype_ht[mt_genes[col_key_name]])
    mt_genes = mt_genes.filter_cols(hl.is_defined(mt_genes[phenotype_field]))

    row_key_name = list(mt_genes.row_key)[0]

    # Collect per-sample qualifying gene sets
    logger.info("Collecting per-sample qualifying gene data...")
    sample_mt = mt_genes.annotate_cols(
        _qualifying_genes=hl.agg.filter(
            mt_genes._qualifies,
            hl.agg.collect_as_set(mt_genes[row_key_name]),
        )
    )

    fields = [phenotype_field, "_qualifying_genes"]
    if covariate_fields:
        fields.extend(covariate_fields)
    samples_df = sample_mt.cols().select(*fields).to_pandas()

    n_samples = len(samples_df)
    logger.info("  %d samples with phenotype data", n_samples)

    if n_samples == 0:
        logger.warning(
            "No samples with phenotype data after joining phenotype table. "
            "Check that sample IDs match."
        )
        return pd.DataFrame(
            columns=[
                "gene_set_name",
                "n_genes_tested",
                "observed_statistic",
                "empirical_p_value",
                "n_permutations",
                "n_exceeded",
            ]
        )

    # Get all background genes
    all_genes_list = mt_genes.aggregate_rows(hl.agg.collect(mt_genes[row_key_name]))
    gene_to_idx = {g: i for i, g in enumerate(sorted(all_genes_list))}
    n_genes = len(gene_to_idx)
    logger.info("  %d background genes", n_genes)

    # Build qualifying matrix: (n_genes, n_samples), sparse
    logger.info("Building qualifying matrix...")
    qual_mat = lil_matrix((n_genes, n_samples), dtype=np.int8)
    for j, qg_set in enumerate(samples_df["_qualifying_genes"]):
        if qg_set is not None:
            for g in qg_set:
                if g in gene_to_idx:
                    qual_mat[gene_to_idx[g], j] = 1
    qual_mat = csr_matrix(qual_mat)

    # Length-matched bins — compute from raw (unnormalized) qualifying counts
    # so bins reflect actual gene variant density, not normalized burden.
    if length_matched:
        gene_qual_counts = np.asarray(qual_mat.sum(axis=1)).flatten()
        from hvantk.algorithms.enrichex.constants import DEFAULT_N_LENGTH_BINS

        n_bins = min(DEFAULT_N_LENGTH_BINS, n_genes // 5)
        if n_bins < 2:
            length_matched = False
            logger.warning("Too few genes for length-matched sampling, using uniform")
        else:
            bin_edges = np.percentile(gene_qual_counts, np.linspace(0, 100, n_bins + 1))
            gene_bins = np.digitize(gene_qual_counts, bin_edges[1:-1])
            bin_to_indices: Dict[int, List[int]] = {}
            for idx in range(n_genes):
                b = int(gene_bins[idx])
                if b not in bin_to_indices:
                    bin_to_indices[b] = []
                bin_to_indices[b].append(idx)

    # Gene-length normalization: weight qualifying entries by 1/gene_length
    if normalize_by_length:
        sorted_genes = sorted(gene_to_idx.keys(), key=lambda g: gene_to_idx[g])
        if gene_lengths is not None:
            logger.info("Normalizing qualifying matrix by provided CDS lengths")
            norm_arr = np.array(
                [gene_lengths.get(g, 1000.0) / 1000.0 for g in sorted_genes],
                dtype=np.float64,
            )
        else:
            logger.info("Normalizing qualifying matrix by variant site count proxy")
            _cohort_rows = cohort_mt.rows()
            gene_site_map = (
                _cohort_rows.group_by(gene=_cohort_rows[gene_field])
                .aggregate(_n_sites=hl.agg.count())
                .to_pandas()
            )
            site_dict = dict(zip(gene_site_map["gene"], gene_site_map["_n_sites"]))
            norm_arr = np.array(
                [float(site_dict.get(g, 1.0)) for g in sorted_genes],
                dtype=np.float64,
            )
        norm_arr = np.maximum(norm_arr, 1e-10)  # avoid division by zero
        from scipy.sparse import diags

        qual_mat = diags(1.0 / norm_arr) @ qual_mat
        qual_mat = csr_matrix(qual_mat)

    # Prepare phenotype and covariate arrays
    y = np.array(samples_df[phenotype_field], dtype=np.float64)
    intercept = np.ones(n_samples, dtype=np.float64)
    if covariate_fields:
        cov_cols = [np.array(samples_df[c], dtype=np.float64) for c in covariate_fields]
        cov_matrix = np.column_stack([intercept] + cov_cols)
    else:
        cov_matrix = intercept.reshape(-1, 1)

    rng = np.random.default_rng(seed)

    # Step 3: Run permutation test for each gene set
    results = []
    # Compute progress milestones (every 25% of permutations)
    _progress_pcts = [25, 50, 75, 100]
    _progress_thresholds = [int(n_permutations * p / 100) for p in _progress_pcts]

    t_perm_total = time.time()

    for gs_name, gs_genes in gene_sets.items():
        logger.info("Testing gene set: %s (%d genes)", gs_name, len(gs_genes))

        gs_indices = np.array([gene_to_idx[g] for g in gs_genes if g in gene_to_idx])
        if len(gs_indices) == 0:
            logger.warning(
                "  No genes from '%s' found in background, skipping", gs_name
            )
            continue

        n_gs_genes = len(gs_indices)
        logger.info("  %d/%d genes found in MT", n_gs_genes, len(gs_genes))

        # Observed burden
        gs_mask = np.zeros(n_genes, dtype=np.float64)
        gs_mask[gs_indices] = 1.0
        observed_burden = np.asarray(gs_mask @ qual_mat).flatten()

        X_obs = np.column_stack([observed_burden, cov_matrix])
        observed_stat = _compute_test_statistic(X_obs, y, phenotype_type)

        if np.isnan(observed_stat):
            logger.warning("  Regression failed for '%s', skipping", gs_name)
            continue

        # Permutations
        n_exceeded = 0
        t_gs_perm = time.time()
        _next_milestone = 0
        for perm_i in range(n_permutations):
            if length_matched:
                perm_indices = _sample_length_matched(
                    gs_indices, gene_bins, bin_to_indices, rng
                )
            else:
                perm_indices = rng.choice(n_genes, size=n_gs_genes, replace=False)

            perm_mask = np.zeros(n_genes, dtype=np.float64)
            perm_mask[perm_indices] = 1.0
            perm_burden = np.asarray(perm_mask @ qual_mat).flatten()

            X_perm = np.column_stack([perm_burden, cov_matrix])
            perm_stat = _compute_test_statistic(X_perm, y, phenotype_type)

            if not np.isnan(perm_stat) and abs(perm_stat) >= abs(observed_stat):
                n_exceeded += 1

            # Log progress at 25% milestones
            if (
                _next_milestone < len(_progress_thresholds)
                and (perm_i + 1) >= _progress_thresholds[_next_milestone]
            ):
                logger.info(
                    "    %d%% permutations complete (%d/%d)",
                    _progress_pcts[_next_milestone],
                    perm_i + 1,
                    n_permutations,
                )
                _next_milestone += 1

        gs_perm_elapsed = time.time() - t_gs_perm

        # +1 correction for empirical p-value
        empirical_p = (n_exceeded + 1) / (n_permutations + 1)

        results.append(
            {
                "gene_set_name": gs_name,
                "n_genes_tested": n_gs_genes,
                "observed_statistic": observed_stat,
                "empirical_p_value": empirical_p,
                "n_permutations": n_permutations,
                "n_exceeded": n_exceeded,
            }
        )

        logger.info(
            "  stat=%.4f, empirical_p=%.4e (%.1fs)",
            observed_stat,
            empirical_p,
            gs_perm_elapsed,
        )

    total_perm_elapsed = time.time() - t_perm_total
    logger.info("\n" + "=" * 60)
    logger.info(
        "COMPETITIVE BURDEN TEST COMPLETE (%.1fs total permutation time)",
        total_perm_elapsed,
    )
    logger.info("=" * 60)

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Phase P: artifact-typed wrappers
# ---------------------------------------------------------------------------


@algorithm(
    name="burden_analysis_artifact",
    backends=[Backend.HAIL],
    inputs={
        "cohort": "ExpressionMatrix",
        "gene_sets": "dict[str, GeneSet]",
        "phenotype": "AnnotationTable",
    },
    outputs={"result": "BurdenAnalysisResult"},
    required_backend="hail",
)
def run_burden_analysis_artifact(
    cohort: "ExpressionMatrix",
    gene_sets: "Dict[str, GeneSet]",
    phenotype: "AnnotationTable",
    **kwargs,
):
    """Phase P artifact-typed wrapper for run_burden_analysis.

    Accepts:
      - cohort: ExpressionMatrix(backend='hail-mt') with genotype data
      - gene_sets: dict of GeneSet collections
      - phenotype: AnnotationTable(backend='hail') with sample phenotypes

    Unwraps via .to_hail_mt() / .to_hail() and delegates to
    run_burden_analysis. The legacy hl.MatrixTable-typed function stays
    as the canonical implementation; this wrapper is the canonical entry
    point for callers consuming the artifact contract.

    All keyword arguments forward to run_burden_analysis (phenotype_field,
    covariate_fields, phenotype_type, gene_field, etc.).
    """
    cohort_mt = cohort.to_hail_mt()
    phenotype_ht = phenotype.to_hail()
    gs_dict = {name: gs.to_list() for name, gs in gene_sets.items()}
    return run_burden_analysis(
        cohort_mt=cohort_mt,
        gene_sets=gs_dict,
        phenotype_ht=phenotype_ht,
        **kwargs,
    )


@algorithm(
    name="stratified_burden_analysis_artifact",
    backends=[Backend.HAIL],
    inputs={
        "cohort": "ExpressionMatrix",
        "gene_sets": "dict[str, GeneSet]",
        "phenotype": "AnnotationTable",
    },
    outputs={"result": "dict[str, BurdenAnalysisResult]"},
    required_backend="hail",
)
def run_stratified_burden_analysis_artifact(
    cohort: "ExpressionMatrix",
    gene_sets: "Dict[str, GeneSet]",
    phenotype: "AnnotationTable",
    variant_classes,
    **kwargs,
):
    """Phase P artifact-typed wrapper for run_stratified_burden_analysis.

    Accepts:
      - cohort: ExpressionMatrix(backend='hail-mt') with genotype data
      - gene_sets: dict of GeneSet collections
      - phenotype: AnnotationTable(backend='hail') with sample phenotypes
      - variant_classes: dict mapping class name to VariantFilter

    Unwraps via .to_hail_mt() / .to_hail() and delegates to
    run_stratified_burden_analysis. The legacy hl.MatrixTable-typed function
    stays as the canonical implementation; this wrapper is the canonical entry
    point for callers consuming the artifact contract.

    All keyword arguments forward to run_stratified_burden_analysis
    (phenotype_field, covariate_fields, phenotype_type, gene_field, etc.).
    """
    cohort_mt = cohort.to_hail_mt()
    phenotype_ht = phenotype.to_hail()
    gs_dict = {name: gs.to_list() for name, gs in gene_sets.items()}
    return run_stratified_burden_analysis(
        cohort_mt=cohort_mt,
        gene_sets=gs_dict,
        phenotype_ht=phenotype_ht,
        variant_classes=variant_classes,
        **kwargs,
    )
