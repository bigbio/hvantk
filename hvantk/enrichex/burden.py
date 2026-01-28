"""
Burden testing module using Hail-native regression.

This module provides functions for computing per-sample burden scores
for gene sets and testing their association with phenotypes using
Hail's native logistic and linear regression functions.

The implementation follows the pattern from scripts/logreg_burden_test.py,
which uses a two-step aggregation:
1. Variants → Genes (per sample)
2. Genes → Gene Sets (per sample)

This keeps all computation in Hail's distributed framework.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

try:  # Optional dependency for burden analysis
    import hail as hl
except ModuleNotFoundError as exc:  # pragma: no cover - depends on env
    hl = None  # type: ignore
    _HAIL_IMPORT_ERROR = exc
else:
    _HAIL_IMPORT_ERROR = None

logger = logging.getLogger(__name__)


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
    min_cadd : Optional[float]
        Minimum CADD score threshold
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
    cadd_field : str
        Field name for CADD score
    consequence_field : str
        Field name for consequence annotation
    """

    max_af: float = 0.01
    min_cadd: Optional[float] = 20.0
    consequences: Optional[List[str]] = None
    pass_only: bool = True
    min_gq: int = 20
    min_dp: int = 10

    # Field name mappings (configurable for different schemas)
    af_field: str = "gnomad_af"
    cadd_field: str = "cadd_phred"
    consequence_field: str = "consequence"

    def to_hail_expr(self, mt: hl.MatrixTable) -> hl.expr.BooleanExpression:
        """Convert filter criteria to Hail boolean expression.

        Parameters
        ----------
        mt : hl.MatrixTable
            MatrixTable to filter

        Returns
        -------
        hl.expr.BooleanExpression
            Combined filter expression
        """
        _require_hail()
        filters = []

        # Allele frequency filter
        if self.max_af is not None:
            if self.af_field in mt.row:
                filters.append(mt[self.af_field] <= self.max_af)
            else:
                logger.warning(f"AF field '{self.af_field}' not found in MT")

        # CADD score filter
        if self.min_cadd is not None:
            if self.cadd_field in mt.row:
                filters.append(mt[self.cadd_field] >= self.min_cadd)
            else:
                logger.warning(f"CADD field '{self.cadd_field}' not found in MT")

        # Consequence filter
        if self.consequences:
            if self.consequence_field in mt.row:
                filters.append(
                    hl.literal(self.consequences).contains(mt[self.consequence_field])
                )
            else:
                logger.warning(
                    f"Consequence field '{self.consequence_field}' not found in MT"
                )

        # PASS filter
        if self.pass_only:
            if "filters" in mt.row:
                filters.append(hl.len(mt.filters) == 0)

        # Genotype quality filters
        if self.min_gq is not None and "GQ" in mt.entry:
            filters.append(mt.GQ >= self.min_gq)

        if self.min_dp is not None and "DP" in mt.entry:
            filters.append(mt.DP >= self.min_dp)

        # Combine all filters
        if filters:
            return hl.all(*filters)
        else:
            return hl.bool(True)


def compute_geneset_burden_mt(
    mt: hl.MatrixTable,
    gene_sets: Dict[str, List[str]],
    gene_field: str = "SYMBOL",
    max_af: float = 0.01,
    min_cadd: Optional[float] = 20.0,
    consequences: Optional[List[str]] = None,
    genotype_aggregation: str = "hets",
) -> hl.MatrixTable:
    """Compute per-sample burden for each gene set, returning a MatrixTable.

    This follows the two-step aggregation pattern from logreg_burden_test.py:
    1. First aggregate variants → genes (per sample)
    2. Then aggregate genes → gene sets (per sample)

    Parameters
    ----------
    mt : hl.MatrixTable
        Annotated variant MatrixTable with genotypes
    gene_sets : Dict[str, List[str]]
        Dictionary mapping gene set names to gene lists
    gene_field : str
        Row field containing gene symbol
    max_af : float
        Maximum allele frequency for qualifying variants
    min_cadd : float, optional
        Minimum CADD score for qualifying variants
    consequences : List[str], optional
        List of qualifying VEP consequences
    genotype_aggregation : str
        One of "hets", "homs", "chets", "homs_chets"

    Returns
    -------
    hl.MatrixTable
        MatrixTable with rows = gene sets, cols = samples, entry = burden (int)

    Examples
    --------
    >>> gene_sets = {
    ...     "microglia": ["TREM2", "CD33", "ABI3"],
    ...     "astrocytes": ["GFAP", "AQP4", "S100B"]
    ... }
    >>> mt_burden = compute_geneset_burden_mt(mt, gene_sets, max_af=0.001)
    >>> mt_burden.count()  # (n_gene_sets, n_samples)
    """
    _require_hail()
    logger.info(f"Computing burden for {len(gene_sets)} gene sets")
    logger.info(f"  Genotype aggregation: {genotype_aggregation}")

    # Validate genotype aggregation method
    valid_methods = ["hets", "homs", "chets", "homs_chets"]
    if genotype_aggregation not in valid_methods:
        raise ValueError(
            f"Invalid genotype_aggregation: {genotype_aggregation}. "
            f"Must be one of: {valid_methods}"
        )

    # Check required field exists
    if gene_field not in mt.row:
        raise ValueError(f"Gene field '{gene_field}' not found in MatrixTable")

    # Apply variant filters
    logger.info("Applying variant filters:")
    filters = []

    if max_af is not None:
        logger.info(f"  Max AF: {max_af}")
        if "gnomad_af" in mt.row:
            filters.append(mt.gnomad_af <= max_af)
        else:
            logger.warning("  gnomad_af field not found, skipping AF filter")

    if min_cadd is not None:
        logger.info(f"  Min CADD: {min_cadd}")
        if "cadd_phred" in mt.row:
            filters.append(mt.cadd_phred >= min_cadd)
        else:
            logger.warning("  cadd_phred field not found, skipping CADD filter")

    if consequences:
        logger.info(f"  Consequences: {consequences}")
        if "consequence" in mt.row:
            filters.append(hl.literal(consequences).contains(mt.consequence))
        else:
            logger.warning("  consequence field not found, skipping consequence filter")

    # Apply filters
    if filters:
        mt = mt.filter_rows(hl.all(*filters))
        n_variants = mt.count_rows()
        logger.info(f"  {n_variants} qualifying variants after filtering")
    else:
        logger.warning("No filters applied")

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
        logger.error("No variants found in gene set genes after filtering!")
        raise ValueError("No variants in gene sets - check gene field and gene IDs")

    # STEP 1: Aggregate variants → genes per sample
    logger.info("Step 1: Aggregating variants to genes per sample...")
    mt_genes = mt.group_rows_by(mt[gene_field]).aggregate(
        hets=hl.agg.count_where(mt.GT.is_het()),
        homs=hl.agg.count_where(mt.GT.is_hom_var()),
        chets=hl.agg.count_where(mt.GT.is_het()) >= 2,
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
    if genotype_aggregation == "hets":
        agg_expr = hl.int(hl.agg.sum(hl.if_else(mt_genes.hets > 0, 1, 0)))
    elif genotype_aggregation == "homs":
        agg_expr = hl.int(hl.agg.sum(hl.if_else(mt_genes.homs > 0, 1, 0)))
    elif genotype_aggregation == "chets":
        agg_expr = hl.int(hl.agg.sum(hl.if_else(mt_genes.chets, 1, 0)))
    elif genotype_aggregation == "homs_chets":
        agg_expr = hl.int(
            hl.agg.sum(hl.if_else(mt_genes.chets | (mt_genes.homs > 0), 1, 0))
        )

    mt_burden = mt_genes.group_rows_by(gene_set_name=mt_genes.gene_set_ids).aggregate(
        burden=agg_expr
    )

    n_gene_sets_final = mt_burden.count_rows()
    n_samples = mt_burden.count_cols()
    logger.info(
        f"Burden matrix created: {n_gene_sets_final} gene sets × {n_samples} samples"
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
    from scripts/logreg_burden_test.py and hvantk/utils/stats.py.

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
    logger.info("Running Hail-native logistic regression")
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

    logger.info(f"Regression complete: {result.count()} gene sets tested")

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
    logger.info("Running Hail-native linear regression")
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

    logger.info(f"Regression complete: {result.count()} gene sets tested")

    return result


def run_burden_analysis(
    cohort_mt: hl.MatrixTable,
    gene_sets: Dict[str, List[str]],
    phenotype_ht: hl.Table,
    phenotype_field: str = "is_case",
    covariate_fields: Optional[List[str]] = None,
    phenotype_type: str = "binary",
    gene_field: str = "SYMBOL",
    max_af: float = 0.01,
    min_cadd: Optional[float] = 20.0,
    consequences: Optional[List[str]] = None,
    genotype_aggregation: str = "hets",
) -> hl.Table:
    """Complete burden analysis pipeline.

    This is the main entry point that orchestrates:
    1. Burden computation per gene set
    2. Phenotype/covariate annotation
    3. Regression testing

    Parameters
    ----------
    cohort_mt : hl.MatrixTable
        Annotated cohort MatrixTable
    gene_sets : Dict[str, List[str]]
        Gene sets to test
    phenotype_ht : hl.Table
        Table with sample phenotypes, keyed by sample ID
    phenotype_field : str
        Field containing phenotype
    covariate_fields : List[str], optional
        Fields to use as covariates
    phenotype_type : str
        "binary" or "continuous"
    gene_field : str
        Row field for gene symbol in cohort_mt
    max_af : float
        Maximum allele frequency filter
    min_cadd : float, optional
        Minimum CADD score filter
    consequences : List[str], optional
        Qualifying consequences
    genotype_aggregation : str
        Genotype aggregation method

    Returns
    -------
    hl.Table
        Results with p-values, odds ratios, etc. for each gene set

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
    mt_burden = compute_geneset_burden_mt(
        mt=cohort_mt,
        gene_sets=gene_sets,
        gene_field=gene_field,
        max_af=max_af,
        min_cadd=min_cadd,
        consequences=consequences,
        genotype_aggregation=genotype_aggregation,
    )

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
        raise ValueError("No samples with phenotype data after join!")

    # Step 3: Run regression
    logger.info(f"\n[3/3] Running {phenotype_type} regression...")

    if phenotype_type == "binary":
        result = logistic_burden_test(
            mt_burden, phenotype_field=phenotype_field, covariates=covariate_fields
        )
    else:
        result = linear_burden_test(
            mt_burden, phenotype_field=phenotype_field, covariates=covariate_fields
        )

    logger.info("\n" + "=" * 60)
    logger.info("BURDEN ANALYSIS COMPLETE")
    logger.info("=" * 60)

    return result
