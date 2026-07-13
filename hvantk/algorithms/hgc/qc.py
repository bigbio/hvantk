"""
Quality Control (QC) module for HGC.

This module provides comprehensive quality control functionality for genomic variant data,
including sample-level and variant-level QC metrics computation, filtering, and table export
based on Hail's QC methods.

The module hands back QC *artifacts* (Hail tables / pandas DataFrames); it does not own
plotting. ``QCMetrics.get_sample_metrics_df()`` / ``get_variant_metrics_df()`` are the stable
entry points for building your own figures (e.g. with matplotlib). The one exception is
``QCMetrics.generate_html_report()``, a canned static HTML report for operational triage.

Functions:
    - compute_sample_qc: Compute sample-level QC metrics
    - compute_variant_qc: Compute variant-level QC metrics
    - compute_full_qc: Compute both sample and variant QC metrics
    - extract_qc_metrics: Extract QC metrics as pandas DataFrames
    - filter_samples_by_qc: Filter samples based on QC thresholds
    - filter_variants_by_qc: Filter variants based on QC thresholds
    - get_qc_summary_stats: Get summary statistics for QC metrics
    - save_qc_metrics: Save QC metrics/tables to disk
"""

import logging
import hail as hl
import pandas as pd
import numpy as np
from typing import Optional, Dict, Tuple, List, Union
from pathlib import Path
from hvantk.core.models.backends import algorithm, Backend

logger = logging.getLogger(__name__)


class QCMetrics:
    """Container class for QC metrics and metadata."""

    def __init__(
        self,
        mt: hl.MatrixTable,
        sample_qc: Optional[hl.Table] = None,
        variant_qc: Optional[hl.Table] = None,
    ):
        self.mt = mt
        self.sample_qc = sample_qc
        self.variant_qc = variant_qc
        self._sample_metrics_df = None
        self._variant_metrics_df = None

    @property
    def has_sample_qc(self) -> bool:
        """Check if sample QC metrics are available."""
        return self.sample_qc is not None

    @property
    def has_variant_qc(self) -> bool:
        """Check if variant QC metrics are available."""
        return self.variant_qc is not None

    def get_sample_metrics_df(self) -> pd.DataFrame:
        """Get sample QC metrics as pandas DataFrame."""
        if self._sample_metrics_df is None and self.has_sample_qc:
            self._sample_metrics_df = self.sample_qc.to_pandas()
        return self._sample_metrics_df

    def get_variant_metrics_df(self) -> pd.DataFrame:
        """Get variant QC metrics as pandas DataFrame."""
        if self._variant_metrics_df is None and self.has_variant_qc:
            self._variant_metrics_df = self.variant_qc.to_pandas()
        return self._variant_metrics_df

    def generate_html_report(self, output_path: Union[str, Path], **kwargs):
        """
        Generate comprehensive HTML QC report.

        Parameters
        ----------
        output_path : str or Path
            Output file path for the HTML report
        **kwargs
            Additional arguments passed to generate_qc_report

        Returns
        -------
        Path
            Path to the generated HTML report

        Example
        -------
        >>> qc_results = compute_full_qc(mt)
        >>> report_path = qc_results.generate_html_report('qc_report.html')
        """
        from hvantk.algorithms.visualization.qc_report import generate_qc_report

        return generate_qc_report(self, output_path, **kwargs)


def _prepare_qc_gt(
    mt: hl.MatrixTable,
    call_field: str = "GT",
    prefer_lpgT: bool = True,
    tmp_field: str = "__qc_gt",
) -> hl.MatrixTable:
    """
    Create a temporary entry field with genotypes aligned to the row alleles and sanitized.

    - Prefer LPGT on split datasets when available (local, biallelic indices).
    - Otherwise use the provided call_field (default GT).
    - Set genotype to missing if any allele index is out of bounds for the current row alleles.
    """
    # Choose base genotype expression
    base_field = call_field
    if prefer_lpgT and ("was_split" in mt.row) and ("LPGT" in mt.entry):
        base_field = "LPGT"

    if base_field not in mt.entry:
        raise ValueError(f"Call field '{base_field}' not found in MatrixTable entries")

    base_gt = mt[base_field]

    # Sanitize: set to missing if any allele index >= len(alleles)
    # Protect against missing genotypes
    invalid_gt = hl.is_defined(base_gt) & hl.any(
        lambda i: base_gt[i] >= hl.len(mt.alleles), hl.range(0, base_gt.ploidy)
    )
    qc_gt = hl.if_else(invalid_gt, hl.missing(hl.tcall), base_gt)

    return mt.annotate_entries(**{tmp_field: qc_gt})


def _with_temp_gt(
    mt: hl.MatrixTable, tmp_field: str = "__qc_gt", backup_field: str = "__orig_GT"
) -> hl.MatrixTable:
    """
    Return a MatrixTable where entry field GT is temporarily set to the sanitized
    genotype stored in tmp_field. If GT exists originally, back it up in backup_field.
    """
    if "GT" in mt.entry:
        return mt.annotate_entries(**{backup_field: mt.GT, "GT": mt[tmp_field]})
    else:
        return mt.annotate_entries(GT=mt[tmp_field])


def _restore_gt(mt: hl.MatrixTable, backup_field: str = "__orig_GT") -> hl.MatrixTable:
    """
    Restore original GT from backup_field if present; otherwise drop GT that was added temporarily.
    Always drop the backup_field if present.
    """
    if backup_field in mt.entry:
        mt = mt.annotate_entries(GT=mt[backup_field])
        mt = mt.drop(backup_field)
    else:
        # GT did not exist originally; keep the sanitized GT as current GT
        pass
    return mt


def compute_sample_qc(
    mt: hl.MatrixTable, name: str = "sample_qc", call_field: str = "GT"
) -> hl.MatrixTable:
    """
    Compute sample-level quality control metrics.

    Computes comprehensive sample QC metrics including:
    - Number of called genotypes
    - Number of heterozygous calls
    - Number of homozygous reference calls
    - Number of homozygous alternate calls
    - Call rate
    - Mean genotype quality (GQ)
    - Mean depth (DP)
    - Transition/transversion ratio

    This function properly handles both:
    - Split multi-allelic variants (after split_multi_hts)
    - Non-split multi-allelic variants
    - Biallelic-only datasets

    Args:
        mt: Input MatrixTable
        name: Name for the sample QC annotation (default: 'sample_qc')
        call_field: Name of the call field to use (default: 'GT')

    Returns:
        MatrixTable with sample QC metrics added to column annotations

    Example:
        >>> mt = hl.read_matrix_table('cohort.mt')
        >>> mt_qc = compute_sample_qc(mt)
        >>> mt_qc.col.sample_qc.describe()
    """
    logger.info("Computing sample-level QC metrics")

    try:
        if call_field not in mt.entry and not (
            ("was_split" in mt.row) and ("LPGT" in mt.entry)
        ):
            raise ValueError(
                f"Call field '{call_field}' not found in MatrixTable entries"
            )

        if "variant_ac" in mt.row:
            logger.debug(
                "Dropping existing variant_ac annotation to avoid indexing issues"
            )
            mt = mt.drop("variant_ac")

        tmp_field = "__qc_gt"
        mt_qc = _prepare_qc_gt(
            mt, call_field=call_field, prefer_lpgT=True, tmp_field=tmp_field
        )

        # Preflight: count invalid GT if any slipped through
        invalid_count = mt_qc.aggregate_entries(
            hl.agg.count_where(
                hl.is_defined(mt_qc[tmp_field])
                & hl.any(
                    lambda i: mt_qc[tmp_field][i] >= hl.len(mt_qc.alleles),
                    hl.range(0, mt_qc[tmp_field].ploidy),
                )
            )
        )
        if invalid_count:
            logger.warning(
                f"Preflight: found {invalid_count} invalid temp genotypes before sample_qc; they will be set missing"
            )

        # Ensure compatibility: temporarily set GT to the sanitized tmp field
        mt_for_qc = _with_temp_gt(mt_qc, tmp_field=tmp_field, backup_field="__orig_GT")

        # Call sample_qc without relying on call_field support
        mt_with_qc = hl.sample_qc(mt_for_qc, name=name)

        # Restore original GT if it existed and drop temps
        mt_with_qc = _restore_gt(mt_with_qc, backup_field="__orig_GT")
        mt_with_qc = mt_with_qc.drop(tmp_field)

        if "GQ" in mt.entry:
            mt_with_qc = mt_with_qc.annotate_cols(
                **{f"{name}_mean_gq": hl.agg.mean(mt_with_qc.GQ)}
            )
        if "DP" in mt.entry:
            mt_with_qc = mt_with_qc.annotate_cols(
                **{f"{name}_mean_dp": hl.agg.mean(mt_with_qc.DP)}
            )

        logger.info(
            f"Successfully computed sample QC metrics for {mt_with_qc.count_cols()} samples"
        )
        return mt_with_qc

    except Exception as e:
        logger.error(f"Failed to compute sample QC: {e}")
        raise


def compute_variant_qc(
    mt: hl.MatrixTable, name: str = "variant_qc", call_field: str = "GT"
) -> hl.MatrixTable:
    """
    Compute variant-level quality control metrics.

    Computes comprehensive variant QC metrics including:
    - Allele count and frequency
    - Number of called samples
    - Call rate
    - Hardy-Weinberg equilibrium test
    - Heterozygosity metrics
    - Mean genotype quality and depth

    This function properly handles both:
    - Split multi-allelic variants (after split_multi_hts)
    - Non-split multi-allelic variants
    - Biallelic-only datasets

    Args:
        mt: Input MatrixTable
        name: Name for the variant QC annotation (default: 'variant_qc')
        call_field: Name of the call field to use (default: 'GT')

    Returns:
        MatrixTable with variant QC metrics added to row annotations

    Example:
        >>> mt = hl.read_matrix_table('cohort.mt')
        >>> mt_qc = compute_variant_qc(mt)
        >>> mt_qc.row.variant_qc.describe()
    """
    logger.info("Computing variant-level QC metrics")

    try:
        if call_field not in mt.entry and not (
            ("was_split" in mt.row) and ("LPGT" in mt.entry)
        ):
            raise ValueError(
                f"Call field '{call_field}' not found in MatrixTable entries"
            )

        if "variant_ac" in mt.row:
            logger.debug("Dropping existing variant_ac annotation")
            mt = mt.drop("variant_ac")

        tmp_field = "__qc_gt"
        mt_qc = _prepare_qc_gt(
            mt, call_field=call_field, prefer_lpgT=True, tmp_field=tmp_field
        )

        # Preflight logging as above
        invalid_count = mt_qc.aggregate_entries(
            hl.agg.count_where(
                hl.is_defined(mt_qc[tmp_field])
                & hl.any(
                    lambda i: mt_qc[tmp_field][i] >= hl.len(mt_qc.alleles),
                    hl.range(0, mt_qc[tmp_field].ploidy),
                )
            )
        )
        if invalid_count:
            logger.warning(
                f"Preflight: found {invalid_count} invalid temp genotypes before variant_qc; they will be set missing"
            )

        mt_for_qc = _with_temp_gt(mt_qc, tmp_field=tmp_field, backup_field="__orig_GT")
        mt_with_qc = hl.variant_qc(mt_for_qc, name=name)
        mt_with_qc = _restore_gt(mt_with_qc, backup_field="__orig_GT")
        mt_with_qc = mt_with_qc.drop(tmp_field)

        if "GQ" in mt.entry:
            mt_with_qc = mt_with_qc.annotate_rows(
                **{f"{name}_mean_gq": hl.agg.mean(mt_with_qc.GQ)}
            )
        if "DP" in mt.entry:
            mt_with_qc = mt_with_qc.annotate_rows(
                **{f"{name}_mean_dp": hl.agg.mean(mt_with_qc.DP)}
            )

        logger.info(
            f"Successfully computed variant QC metrics for {mt_with_qc.count_rows()} variants"
        )
        return mt_with_qc

    except Exception as e:
        logger.error(f"Failed to compute variant QC: {e}")
        raise


@algorithm(name="compute_full_qc", backends=[Backend.HAIL])
def compute_full_qc(
    mt: hl.MatrixTable,
    sample_qc_name: str = "sample_qc",
    variant_qc_name: str = "variant_qc",
    call_field: str = "GT",
) -> QCMetrics:
    """
    Compute comprehensive QC metrics for both samples and variants.

    Args:
        mt: Input MatrixTable
        sample_qc_name: Name for sample QC annotation
        variant_qc_name: Name for variant QC annotation
        call_field: Name of the call field to use

    Returns:
        QCMetrics object containing the MatrixTable with QC metrics
        and extracted QC tables

    Example:
        >>> mt = hl.read_matrix_table('cohort.mt')
        >>> qc_results = compute_full_qc(mt)
        >>> sample_df = qc_results.get_sample_metrics_df()
        >>> variant_df = qc_results.get_variant_metrics_df()

    Notes:
        This algorithm operates on raw `hl.MatrixTable` / `hl.VariantDataset` instances
        (genotype data). ExpressionMatrix's hail-mt backend isn't available yet (Phase J).
    """
    logger.info("Computing comprehensive QC metrics for samples and variants")

    try:
        # Compute sample QC
        mt_with_sample_qc = compute_sample_qc(
            mt, name=sample_qc_name, call_field=call_field
        )

        # Compute variant QC
        mt_with_full_qc = compute_variant_qc(
            mt_with_sample_qc, name=variant_qc_name, call_field=call_field
        )

        # Extract QC tables
        sample_qc_table = mt_with_full_qc.cols().select(sample_qc_name)
        variant_qc_table = mt_with_full_qc.rows().select(variant_qc_name)

        logger.info("Successfully computed full QC metrics")
        return QCMetrics(mt_with_full_qc, sample_qc_table, variant_qc_table)

    except Exception as e:
        logger.error(f"Failed to compute full QC: {e}")
        raise


def extract_qc_metrics(
    mt: hl.MatrixTable,
    sample_qc_name: str = "sample_qc",
    variant_qc_name: str = "variant_qc",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Extract QC metrics as pandas DataFrames.

    Args:
        mt: MatrixTable with QC annotations
        sample_qc_name: Name of sample QC annotation
        variant_qc_name: Name of variant QC annotation

    Returns:
        Tuple of (sample_qc_df, variant_qc_df)

    Example:
        >>> mt_qc = compute_full_qc(mt)
        >>> sample_df, variant_df = extract_qc_metrics(mt_qc.mt)
    """
    logger.info("Extracting QC metrics as pandas DataFrames")

    try:
        sample_qc_df = None
        variant_qc_df = None

        # Extract sample QC if available
        if sample_qc_name in mt.col:
            sample_qc_table = mt.cols().select(sample_qc_name)
            sample_qc_df = sample_qc_table.to_pandas()
            logger.info(f"Extracted sample QC metrics for {len(sample_qc_df)} samples")

        # Extract variant QC if available
        if variant_qc_name in mt.row:
            variant_qc_table = mt.rows().select(variant_qc_name)
            variant_qc_df = variant_qc_table.to_pandas()
            logger.info(
                f"Extracted variant QC metrics for {len(variant_qc_df)} variants"
            )

        return sample_qc_df, variant_qc_df

    except Exception as e:
        logger.error(f"Failed to extract QC metrics: {e}")
        raise


def filter_samples_by_qc(
    mt: hl.MatrixTable,
    min_call_rate: float = 0.85,
    min_mean_dp: Optional[float] = None,
    max_mean_dp: Optional[float] = None,
    min_mean_gq: Optional[float] = None,
    sample_qc_name: str = "sample_qc",
) -> hl.MatrixTable:
    """
    Filter samples based on QC thresholds.

    Args:
        mt: MatrixTable with sample QC annotations
        min_call_rate: Minimum sample call rate (default: 0.85)
        min_mean_dp: Minimum mean depth (optional)
        max_mean_dp: Maximum mean depth (optional)
        min_mean_gq: Minimum mean genotype quality (optional)
        sample_qc_name: Name of sample QC annotation

    Returns:
        Filtered MatrixTable

    Example:
        >>> mt_qc = compute_sample_qc(mt)
        >>> mt_filtered = filter_samples_by_qc(mt_qc, min_call_rate=0.9, min_mean_dp=10)
    """
    logger.info("Filtering samples based on QC thresholds")

    try:
        if sample_qc_name not in mt.col:
            raise ValueError(f"Sample QC annotation '{sample_qc_name}' not found")

        # Build filter expression
        filters = []

        # Call rate filter
        filters.append(mt[sample_qc_name].call_rate >= min_call_rate)

        # Depth filters
        if min_mean_dp is not None and f"{sample_qc_name}_mean_dp" in mt.col:
            filters.append(mt[f"{sample_qc_name}_mean_dp"] >= min_mean_dp)

        if max_mean_dp is not None and f"{sample_qc_name}_mean_dp" in mt.col:
            filters.append(mt[f"{sample_qc_name}_mean_dp"] <= max_mean_dp)

        # Genotype quality filter
        if min_mean_gq is not None and f"{sample_qc_name}_mean_gq" in mt.col:
            filters.append(mt[f"{sample_qc_name}_mean_gq"] >= min_mean_gq)

        # Apply combined filter
        filter_expr = hl.all(filters) if len(filters) > 1 else filters[0]
        mt_filtered = mt.filter_cols(filter_expr)

        n_samples_before = mt.count_cols()
        n_samples_after = mt_filtered.count_cols()
        n_removed = n_samples_before - n_samples_after

        logger.info(
            f"Filtered {n_removed} samples ({n_removed/n_samples_before*100:.1f}%), "
            f"{n_samples_after} samples remaining"
        )

        return mt_filtered

    except Exception as e:
        logger.error(f"Failed to filter samples: {e}")
        raise


def filter_variants_by_qc(
    mt: hl.MatrixTable,
    min_call_rate: float = 0.85,
    min_ac: int = 1,
    max_ac: Optional[int] = None,
    min_af: Optional[float] = None,
    max_af: Optional[float] = None,
    hwe_threshold: Optional[float] = 1e-6,
    variant_qc_name: str = "variant_qc",
) -> hl.MatrixTable:
    """
    Filter variants based on QC thresholds.

    Args:
        mt: MatrixTable with variant QC annotations
        min_call_rate: Minimum variant call rate (default: 0.85)
        min_ac: Minimum allele count (default: 1)
        max_ac: Maximum allele count (optional)
        min_af: Minimum allele frequency (optional)
        max_af: Maximum allele frequency (optional)
        hwe_threshold: Hardy-Weinberg equilibrium p-value threshold (optional)
        variant_qc_name: Name of variant QC annotation

    Returns:
        Filtered MatrixTable

    Example:
        >>> mt_qc = compute_variant_qc(mt)
        >>> mt_filtered = filter_variants_by_qc(mt_qc, min_call_rate=0.9, min_ac=5)
    """
    logger.info("Filtering variants based on QC thresholds")

    try:
        if variant_qc_name not in mt.row:
            raise ValueError(f"Variant QC annotation '{variant_qc_name}' not found")

        # Build filter expression
        filters = []

        # Call rate filter
        filters.append(mt[variant_qc_name].call_rate >= min_call_rate)

        # Allele count filters
        filters.append(mt[variant_qc_name].AC[1] >= min_ac)

        if max_ac is not None:
            filters.append(mt[variant_qc_name].AC[1] <= max_ac)

        # Allele frequency filters
        if min_af is not None:
            filters.append(mt[variant_qc_name].AF[1] >= min_af)

        if max_af is not None:
            filters.append(mt[variant_qc_name].AF[1] <= max_af)

        # Hardy-Weinberg equilibrium filter
        if hwe_threshold is not None:
            filters.append(mt[variant_qc_name].p_value_hwe >= hwe_threshold)

        # Apply combined filter
        filter_expr = hl.all(filters) if len(filters) > 1 else filters[0]
        mt_filtered = mt.filter_rows(filter_expr)

        n_variants_before = mt.count_rows()
        n_variants_after = mt_filtered.count_rows()
        n_removed = n_variants_before - n_variants_after

        logger.info(
            f"Filtered {n_removed} variants ({n_removed/n_variants_before*100:.1f}%), "
            f"{n_variants_after} variants remaining"
        )

        return mt_filtered

    except Exception as e:
        logger.error(f"Failed to filter variants: {e}")
        raise


def get_qc_summary_stats(
    qc_df: pd.DataFrame, metrics: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Get summary statistics for QC metrics.

    Args:
        qc_df: DataFrame with QC metrics
        metrics: List of specific metrics to summarize (optional)

    Returns:
        DataFrame with summary statistics

    Example:
        >>> sample_df, _ = extract_qc_metrics(mt_qc)
        >>> stats = get_qc_summary_stats(sample_df, ['call_rate', 'n_called'])
    """
    try:
        if metrics is None:
            # Auto-detect numeric columns
            numeric_cols = qc_df.select_dtypes(include=[np.number]).columns.tolist()
        else:
            numeric_cols = [col for col in metrics if col in qc_df.columns]

        if not numeric_cols:
            logger.warning("No numeric QC metrics found for summary statistics")
            return pd.DataFrame()

        summary = qc_df[numeric_cols].describe()
        logger.info(f"Generated summary statistics for {len(numeric_cols)} QC metrics")

        return summary

    except Exception as e:
        logger.error(f"Failed to generate QC summary statistics: {e}")
        raise


def prepare_qc_for_visualization(
    qc_metrics: QCMetrics,
    sample_metrics: Optional[List[str]] = None,
    variant_metrics: Optional[List[str]] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Prepare QC data for visualization by extracting and flattening key metrics.

    Args:
        qc_metrics: QCMetrics object with computed QC data
        sample_metrics: List of sample metrics to prepare (optional)
        variant_metrics: List of variant metrics to prepare (optional)

    Returns:
        Dictionary with 'sample' and 'variant' DataFrames ready for plotting

    Example:
        >>> qc_results = compute_full_qc(mt)
        >>> viz_data = prepare_qc_for_visualization(qc_results)
        >>> sample_data = viz_data['sample']
        >>> variant_data = viz_data['variant']
    """
    logger.info("Preparing QC data for visualization")

    result = {}

    try:
        # Prepare sample QC data
        if qc_metrics.has_sample_qc:
            sample_df = qc_metrics.get_sample_metrics_df()

            # Flatten nested column names for easier plotting access
            flattened_sample = sample_df.copy()
            new_columns = []
            for col in flattened_sample.columns:
                if isinstance(col, str) and "." in col:
                    # Extract the metric name after the last dot
                    new_name = col.split(".")[-1]
                    new_columns.append(new_name)
                else:
                    new_columns.append(col)

            flattened_sample.columns = new_columns

            if sample_metrics is None:
                # Use all available flattened metrics
                result["sample"] = flattened_sample
            else:
                # Filter to requested metrics
                available_sample_metrics = [
                    m for m in sample_metrics if m in flattened_sample.columns
                ]
                result["sample"] = flattened_sample[available_sample_metrics]

            logger.info(
                f"Prepared {len(result['sample'].columns)} sample metrics for visualization"
            )

        # Prepare variant QC data
        if qc_metrics.has_variant_qc:
            variant_df = qc_metrics.get_variant_metrics_df()

            # Flatten nested column names and handle array columns
            flattened_variant = variant_df.copy()
            new_columns = []
            for col in flattened_variant.columns:
                if isinstance(col, str) and "." in col:
                    new_name = col.split(".")[-1]
                    new_columns.append(new_name)
                else:
                    new_columns.append(col)

            flattened_variant.columns = new_columns

            # Handle array columns (AC, AF) by extracting alternate allele values
            if "AC" in flattened_variant.columns:
                try:
                    flattened_variant["AC_alt"] = flattened_variant["AC"].apply(
                        lambda x: (
                            x[1]
                            if isinstance(x, (list, np.ndarray)) and len(x) > 1
                            else x
                        )
                    )
                except (TypeError, IndexError):
                    logger.warning(
                        "Could not extract alternate allele count from AC column"
                    )

            if "AF" in flattened_variant.columns:
                try:
                    flattened_variant["AF_alt"] = flattened_variant["AF"].apply(
                        lambda x: (
                            x[1]
                            if isinstance(x, (list, np.ndarray)) and len(x) > 1
                            else x
                        )
                    )
                except (TypeError, IndexError):
                    logger.warning(
                        "Could not extract alternate allele frequency from AF column"
                    )

            if variant_metrics is None:
                # Use all available flattened metrics
                result["variant"] = flattened_variant
            else:
                # Filter to requested metrics
                available_variant_metrics = [
                    m for m in variant_metrics if m in flattened_variant.columns
                ]
                result["variant"] = flattened_variant[available_variant_metrics]

            logger.info(
                f"Prepared {len(result['variant'].columns)} variant metrics for visualization"
            )

        return result

    except Exception as e:
        logger.error(f"Failed to prepare QC data for visualization: {e}")
        raise


def save_qc_metrics(
    qc_metrics: QCMetrics, output_dir: Union[str, Path], prefix: str = "qc_metrics"
) -> Dict[str, str]:
    """
    Save QC metrics to files.

    Args:
        qc_metrics: QCMetrics object
        output_dir: Output directory
        prefix: File prefix for output files

    Returns:
        Dictionary with paths to saved files

    Example:
        >>> qc_results = compute_full_qc(mt)
        >>> saved_files = save_qc_metrics(qc_results, 'qc_output/')
    """
    logger.info("Saving QC metrics to files")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    saved_files = {}

    try:
        # Save sample QC
        if qc_metrics.has_sample_qc:
            sample_df = qc_metrics.get_sample_metrics_df()
            sample_path = output_dir / f"{prefix}_sample_qc.csv"
            sample_df.to_csv(sample_path, index=False)
            saved_files["sample_qc"] = str(sample_path)
            logger.info(f"Saved sample QC metrics to {sample_path}")

        # Save variant QC
        if qc_metrics.has_variant_qc:
            variant_df = qc_metrics.get_variant_metrics_df()
            variant_path = output_dir / f"{prefix}_variant_qc.csv"
            variant_df.to_csv(variant_path, index=False)
            saved_files["variant_qc"] = str(variant_path)
            logger.info(f"Saved variant QC metrics to {variant_path}")

        # Save MatrixTable if needed
        mt_path = output_dir / f"{prefix}_with_qc.mt"
        qc_metrics.mt.write(str(mt_path), overwrite=True)
        saved_files["matrix_table"] = str(mt_path)
        logger.info(f"Saved MatrixTable with QC annotations to {mt_path}")

        return saved_files

    except Exception as e:
        logger.error(f"Failed to save QC metrics: {e}")
        raise
