"""
Quality Control (QC) module for HGC.

This module provides comprehensive quality control functionality for genomic variant data,
including sample-level and variant-level QC metrics computation, filtering, and visualization
preparation based on Hail's QC methods.

Functions:
    - compute_sample_qc: Compute sample-level QC metrics
    - compute_variant_qc: Compute variant-level QC metrics
    - compute_full_qc: Compute both sample and variant QC metrics
    - extract_qc_metrics: Extract QC metrics as pandas DataFrames
    - filter_samples_by_qc: Filter samples based on QC thresholds
    - filter_variants_by_qc: Filter variants based on QC thresholds
    - get_qc_summary_stats: Get summary statistics for QC metrics
    - prepare_qc_for_visualization: Prepare QC data for plotting
"""

import logging
import hail as hl
import pandas as pd
import numpy as np
from typing import Optional, Dict, Tuple, List, Union
from pathlib import Path

logger = logging.getLogger(__name__)


class QCMetrics:
    """Container class for QC metrics and metadata."""

    def __init__(self, mt: hl.MatrixTable, sample_qc: Optional[hl.Table] = None,
                 variant_qc: Optional[hl.Table] = None):
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

    def plot_sample_overview(self, **kwargs):
        """
        Plot sample QC overview using hvantk visualization functions.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_sample_qc_overview

        Returns
        -------
        matplotlib.figure.Figure
            The created figure
        """
        if not self.has_sample_qc:
            raise ValueError("No sample QC data available for plotting")

        from hvantk.visualization.qc_plots import plot_sample_qc_overview
        return plot_sample_qc_overview(self.get_sample_metrics_df(), **kwargs)

    def plot_variant_overview(self, **kwargs):
        """
        Plot variant QC overview using hvantk visualization functions.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_variant_qc_overview

        Returns
        -------
        matplotlib.figure.Figure
            The created figure
        """
        if not self.has_variant_qc:
            raise ValueError("No variant QC data available for plotting")

        from hvantk.visualization.qc_plots import plot_variant_qc_overview
        return plot_variant_qc_overview(self.get_variant_metrics_df(), **kwargs)

    def plot_dashboard(self, **kwargs):
        """
        Plot comprehensive QC summary dashboard.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_qc_summary_dashboard

        Returns
        -------
        matplotlib.figure.Figure
            The created figure
        """
        from hvantk.visualization.qc_plots import plot_qc_summary_dashboard
        return plot_qc_summary_dashboard(self, **kwargs)

    def plot_sample_call_rates(self, **kwargs):
        """
        Plot sample call rate distribution.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_sample_call_rate_distribution

        Returns
        -------
        matplotlib.figure.Figure
            The created figure
        """
        if not self.has_sample_qc:
            raise ValueError("No sample QC data available for plotting")

        from hvantk.visualization.qc_plots import plot_sample_call_rate_distribution
        return plot_sample_call_rate_distribution(self.get_sample_metrics_df(), **kwargs)

    def plot_sample_titv(self, **kwargs):
        """
        Plot sample Ti/Tv ratio distribution.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_sample_titv_distribution

        Returns
        -------
        matplotlib.figure.Figure
            The created figure
        """
        if not self.has_sample_qc:
            raise ValueError("No sample QC data available for plotting")

        from hvantk.visualization.qc_plots import plot_sample_titv_distribution
        return plot_sample_titv_distribution(self.get_sample_metrics_df(), **kwargs)

    def plot_variant_call_rates(self, **kwargs):
        """
        Plot variant call rate distribution.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_variant_call_rate_distribution

        Returns
        -------
        matplotlib.figure.Figure
            The created figure
        """
        if not self.has_variant_qc:
            raise ValueError("No variant QC data available for plotting")

        from hvantk.visualization.qc_plots import plot_variant_call_rate_distribution
        return plot_variant_call_rate_distribution(self.get_variant_metrics_df(), **kwargs)

    def plot_allele_frequencies(self, **kwargs):
        """
        Plot allele frequency spectrum.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_allele_frequency_spectrum

        Returns
        -------
        matplotlib.figure.Figure
            The created figure
        """
        if not self.has_variant_qc:
            raise ValueError("No variant QC data available for plotting")

        from hvantk.visualization.qc_plots import plot_allele_frequency_spectrum
        return plot_allele_frequency_spectrum(self.get_variant_metrics_df(), **kwargs)

    def plot_hwe_pvalues(self, **kwargs):
        """
        Plot Hardy-Weinberg equilibrium p-values.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_hwe_pvalues

        Returns
        -------
        matplotlib.figure.Figure
            The created figure
        """
        if not self.has_variant_qc:
            raise ValueError("No variant QC data available for plotting")

        from hvantk.visualization.qc_plots import plot_hwe_pvalues
        return plot_hwe_pvalues(self.get_variant_metrics_df(), **kwargs)

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
        from hvantk.visualization.qc_report import generate_qc_report
        return generate_qc_report(self, output_path, **kwargs)

    # Interactive plotting methods
    def plot_interactive_sample_call_rates(self, **kwargs):
        """
        Plot interactive sample call rate distribution.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_interactive_sample_call_rates

        Returns
        -------
        plotly.graph_objects.Figure
            Interactive plotly figure
        """
        if not self.has_sample_qc:
            raise ValueError("No sample QC data available for plotting")

        from hvantk.visualization.interactive_qc import plot_interactive_sample_call_rates
        return plot_interactive_sample_call_rates(self.get_sample_metrics_df(), **kwargs)

    def plot_interactive_sample_titv(self, **kwargs):
        """
        Plot interactive sample Ti/Tv ratio distribution.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_interactive_sample_titv

        Returns
        -------
        plotly.graph_objects.Figure
            Interactive plotly figure
        """
        if not self.has_sample_qc:
            raise ValueError("No sample QC data available for plotting")

        from hvantk.visualization.interactive_qc import plot_interactive_sample_titv
        return plot_interactive_sample_titv(self.get_sample_metrics_df(), **kwargs)

    def plot_interactive_variant_call_rates(self, **kwargs):
        """
        Plot interactive variant call rate distribution.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_interactive_variant_call_rates

        Returns
        -------
        plotly.graph_objects.Figure
            Interactive plotly figure
        """
        if not self.has_variant_qc:
            raise ValueError("No variant QC data available for plotting")

        from hvantk.visualization.interactive_qc import plot_interactive_variant_call_rates
        return plot_interactive_variant_call_rates(self.get_variant_metrics_df(), **kwargs)

    def plot_interactive_allele_frequencies(self, **kwargs):
        """
        Plot interactive allele frequency spectrum.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_interactive_allele_frequencies

        Returns
        -------
        plotly.graph_objects.Figure
            Interactive plotly figure
        """
        if not self.has_variant_qc:
            raise ValueError("No variant QC data available for plotting")

        from hvantk.visualization.interactive_qc import plot_interactive_allele_frequencies
        return plot_interactive_allele_frequencies(self.get_variant_metrics_df(), **kwargs)

    def plot_interactive_hwe_pvalues(self, **kwargs):
        """
        Plot interactive Hardy-Weinberg equilibrium p-values.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_interactive_hwe_pvalues

        Returns
        -------
        plotly.graph_objects.Figure
            Interactive plotly figure
        """
        if not self.has_variant_qc:
            raise ValueError("No variant QC data available for plotting")

        from hvantk.visualization.interactive_qc import plot_interactive_hwe_pvalues
        return plot_interactive_hwe_pvalues(self.get_variant_metrics_df(), **kwargs)

    def plot_interactive_dashboard(self, **kwargs):
        """
        Plot interactive comprehensive QC dashboard.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_interactive_qc_dashboard

        Returns
        -------
        plotly.graph_objects.Figure
            Interactive plotly figure
        """
        from hvantk.visualization.interactive_qc import plot_interactive_qc_dashboard
        return plot_interactive_qc_dashboard(self, **kwargs)

    def plot_interactive_sample_scatter(self, x_metric='call_rate', y_metric='r_ti_tv', **kwargs):
        """
        Plot interactive sample scatter plot.

        Parameters
        ----------
        x_metric : str
            Column name for x-axis metric
        y_metric : str
            Column name for y-axis metric
        **kwargs
            Additional arguments passed to plot_interactive_sample_scatter

        Returns
        -------
        plotly.graph_objects.Figure
            Interactive plotly figure
        """
        if not self.has_sample_qc:
            raise ValueError("No sample QC data available for plotting")

        from hvantk.visualization.interactive_qc import plot_interactive_sample_scatter
        return plot_interactive_sample_scatter(
            self.get_sample_metrics_df(),
            x_metric=x_metric,
            y_metric=y_metric,
            **kwargs
        )




def compute_sample_qc(mt: hl.MatrixTable,
                     name: str = 'sample_qc',
                     call_field: str = 'GT') -> hl.MatrixTable:
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
        # Ensure the call field exists
        if call_field not in mt.entry:
            raise ValueError(f"Call field '{call_field}' not found in MatrixTable entries")

        # Remove any existing variant_ac annotation to avoid Hail's internal
        # array indexing bug with split multi-allelic variants.
        # Hail's sample_qc() will recompute variant_ac internally if needed.
        if "variant_ac" in mt.row:
            logger.debug("Dropping existing variant_ac annotation to avoid indexing issues")
            mt = mt.drop("variant_ac")

        # CRITICAL FIX: For split multi-allelic variants, use LPGT instead of GT
        # GT may still have original allele indices, but LPGT has local (biallelic) indices
        effective_call_field = call_field
        had_gt_originally = 'GT' in mt.entry

        # Detect if variants have been split and LPGT is available
        if 'was_split' in mt.row and 'LPGT' in mt.entry:
            logger.info("Detected split multi-allelic variants - using LPGT for QC to avoid array indexing issues")
            effective_call_field = 'LPGT'

        # If using a non-standard call field, we need to ensure Hail uses it via GT
        if effective_call_field != 'GT':
            logger.info(f"Using '{effective_call_field}' as the genotype field for sample QC")
            # Temporarily replace GT with the effective call field
            mt = mt.annotate_entries(GT=mt[effective_call_field])

        # Compute sample QC using Hail's built-in function
        # This will internally compute variant_ac correctly for the current genotypes
        mt_with_qc = hl.sample_qc(mt, name=name)

        # Restore original GT state: if we used a different field, restore GT from effective_call_field
        # This ensures GT is bound to the new MatrixTable source
        if effective_call_field != 'GT':
            if had_gt_originally:
                # Re-annotate GT from the effective call field (now from new MT source)
                mt_with_qc = mt_with_qc.annotate_entries(GT=mt_with_qc[effective_call_field])
            else:
                # GT didn't exist originally, so drop it
                mt_with_qc = mt_with_qc.drop('GT')

        # Add additional custom metrics if available
        if 'GQ' in mt.entry:
            # Add mean genotype quality
            mt_with_qc = mt_with_qc.annotate_cols(
                **{f"{name}_mean_gq": hl.agg.mean(mt_with_qc.GQ)}
            )

        if 'DP' in mt.entry:
            # Add mean depth
            mt_with_qc = mt_with_qc.annotate_cols(
                **{f"{name}_mean_dp": hl.agg.mean(mt_with_qc.DP)}
            )

        logger.info(f"Successfully computed sample QC metrics for {mt_with_qc.count_cols()} samples")
        return mt_with_qc

    except Exception as e:
        logger.error(f"Failed to compute sample QC: {e}")
        raise


def compute_variant_qc(mt: hl.MatrixTable,
                      name: str = 'variant_qc',
                      call_field: str = 'GT') -> hl.MatrixTable:
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
        # Ensure the call field exists
        if call_field not in mt.entry:
            raise ValueError(f"Call field '{call_field}' not found in MatrixTable entries")

        # Remove any existing variant_ac to let variant_qc compute it fresh
        if "variant_ac" in mt.row:
            logger.debug("Dropping existing variant_ac annotation")
            mt = mt.drop("variant_ac")

        # CRITICAL FIX: For split multi-allelic variants, use LPGT instead of GT
        # GT may still have original allele indices, but LPGT has local (biallelic) indices
        effective_call_field = call_field
        had_gt_originally = 'GT' in mt.entry

        # Detect if variants have been split and LPGT is available
        if 'was_split' in mt.row and 'LPGT' in mt.entry:
            logger.info("Detected split multi-allelic variants - using LPGT for QC to avoid array indexing issues")
            effective_call_field = 'LPGT'

        # If using a non-standard call field, ensure Hail uses it via GT
        if effective_call_field != 'GT':
            logger.info(f"Using '{effective_call_field}' as the genotype field for variant QC")
            # Temporarily replace GT with the effective call field
            mt = mt.annotate_entries(GT=mt[effective_call_field])

        # Compute variant QC using Hail's built-in function
        mt_with_qc = hl.variant_qc(mt, name=name)

        # Restore original GT state: if we used a different field, restore GT from effective_call_field
        # This ensures GT is bound to the new MatrixTable source
        if effective_call_field != 'GT':
            if had_gt_originally:
                # Re-annotate GT from the effective call field (now from new MT source)
                mt_with_qc = mt_with_qc.annotate_entries(GT=mt_with_qc[effective_call_field])
            else:
                # GT didn't exist originally, so drop it
                mt_with_qc = mt_with_qc.drop('GT')

        # Add additional custom metrics if available
        if 'GQ' in mt.entry:
            # Add mean genotype quality per variant
            mt_with_qc = mt_with_qc.annotate_rows(
                **{f"{name}_mean_gq": hl.agg.mean(mt_with_qc.GQ)}
            )

        if 'DP' in mt.entry:
            # Add mean depth per variant
            mt_with_qc = mt_with_qc.annotate_rows(
                **{f"{name}_mean_dp": hl.agg.mean(mt_with_qc.DP)}
            )

        logger.info(f"Successfully computed variant QC metrics for {mt_with_qc.count_rows()} variants")
        return mt_with_qc

    except Exception as e:
        logger.error(f"Failed to compute variant QC: {e}")
        raise


def compute_full_qc(mt: hl.MatrixTable,
                   sample_qc_name: str = 'sample_qc',
                   variant_qc_name: str = 'variant_qc',
                   call_field: str = 'GT') -> QCMetrics:
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
    """
    logger.info("Computing comprehensive QC metrics for samples and variants")

    try:
        # Compute sample QC
        mt_with_sample_qc = compute_sample_qc(mt, name=sample_qc_name, call_field=call_field)

        # Compute variant QC
        mt_with_full_qc = compute_variant_qc(mt_with_sample_qc, name=variant_qc_name, call_field=call_field)

        # Extract QC tables
        sample_qc_table = mt_with_full_qc.cols().select(sample_qc_name)
        variant_qc_table = mt_with_full_qc.rows().select(variant_qc_name)

        logger.info("Successfully computed full QC metrics")
        return QCMetrics(mt_with_full_qc, sample_qc_table, variant_qc_table)

    except Exception as e:
        logger.error(f"Failed to compute full QC: {e}")
        raise


def extract_qc_metrics(mt: hl.MatrixTable,
                      sample_qc_name: str = 'sample_qc',
                      variant_qc_name: str = 'variant_qc') -> Tuple[pd.DataFrame, pd.DataFrame]:
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
            logger.info(f"Extracted variant QC metrics for {len(variant_qc_df)} variants")

        return sample_qc_df, variant_qc_df

    except Exception as e:
        logger.error(f"Failed to extract QC metrics: {e}")
        raise


def filter_samples_by_qc(mt: hl.MatrixTable,
                        min_call_rate: float = 0.85,
                        min_mean_dp: Optional[float] = None,
                        max_mean_dp: Optional[float] = None,
                        min_mean_gq: Optional[float] = None,
                        sample_qc_name: str = 'sample_qc') -> hl.MatrixTable:
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

        logger.info(f"Filtered {n_removed} samples ({n_removed/n_samples_before*100:.1f}%), "
                   f"{n_samples_after} samples remaining")

        return mt_filtered

    except Exception as e:
        logger.error(f"Failed to filter samples: {e}")
        raise


def filter_variants_by_qc(mt: hl.MatrixTable,
                         min_call_rate: float = 0.85,
                         min_ac: int = 1,
                         max_ac: Optional[int] = None,
                         min_af: Optional[float] = None,
                         max_af: Optional[float] = None,
                         hwe_threshold: Optional[float] = 1e-6,
                         variant_qc_name: str = 'variant_qc') -> hl.MatrixTable:
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

        logger.info(f"Filtered {n_removed} variants ({n_removed/n_variants_before*100:.1f}%), "
                   f"{n_variants_after} variants remaining")

        return mt_filtered

    except Exception as e:
        logger.error(f"Failed to filter variants: {e}")
        raise


def get_qc_summary_stats(qc_df: pd.DataFrame,
                        metrics: Optional[List[str]] = None) -> pd.DataFrame:
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


def prepare_qc_for_visualization(qc_metrics: QCMetrics,
                                sample_metrics: Optional[List[str]] = None,
                                variant_metrics: Optional[List[str]] = None) -> Dict[str, pd.DataFrame]:
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
                if isinstance(col, str) and '.' in col:
                    # Extract the metric name after the last dot
                    new_name = col.split('.')[-1]
                    new_columns.append(new_name)
                else:
                    new_columns.append(col)

            flattened_sample.columns = new_columns

            if sample_metrics is None:
                # Use all available flattened metrics
                result['sample'] = flattened_sample
            else:
                # Filter to requested metrics
                available_sample_metrics = [m for m in sample_metrics if m in flattened_sample.columns]
                result['sample'] = flattened_sample[available_sample_metrics]

            logger.info(f"Prepared {len(result['sample'].columns)} sample metrics for visualization")

        # Prepare variant QC data
        if qc_metrics.has_variant_qc:
            variant_df = qc_metrics.get_variant_metrics_df()

            # Flatten nested column names and handle array columns
            flattened_variant = variant_df.copy()
            new_columns = []
            for col in flattened_variant.columns:
                if isinstance(col, str) and '.' in col:
                    new_name = col.split('.')[-1]
                    new_columns.append(new_name)
                else:
                    new_columns.append(col)

            flattened_variant.columns = new_columns

            # Handle array columns (AC, AF) by extracting alternate allele values
            if 'AC' in flattened_variant.columns:
                try:
                    flattened_variant['AC_alt'] = flattened_variant['AC'].apply(
                        lambda x: x[1] if isinstance(x, (list, np.ndarray)) and len(x) > 1 else x
                    )
                except (TypeError, IndexError):
                    logger.warning("Could not extract alternate allele count from AC column")

            if 'AF' in flattened_variant.columns:
                try:
                    flattened_variant['AF_alt'] = flattened_variant['AF'].apply(
                        lambda x: x[1] if isinstance(x, (list, np.ndarray)) and len(x) > 1 else x
                    )
                except (TypeError, IndexError):
                    logger.warning("Could not extract alternate allele frequency from AF column")

            if variant_metrics is None:
                # Use all available flattened metrics
                result['variant'] = flattened_variant
            else:
                # Filter to requested metrics
                available_variant_metrics = [m for m in variant_metrics if m in flattened_variant.columns]
                result['variant'] = flattened_variant[available_variant_metrics]

            logger.info(f"Prepared {len(result['variant'].columns)} variant metrics for visualization")

        return result

    except Exception as e:
        logger.error(f"Failed to prepare QC data for visualization: {e}")
        raise


def save_qc_metrics(qc_metrics: QCMetrics,
                   output_dir: Union[str, Path],
                   prefix: str = 'qc_metrics') -> Dict[str, str]:
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
            saved_files['sample_qc'] = str(sample_path)
            logger.info(f"Saved sample QC metrics to {sample_path}")

        # Save variant QC
        if qc_metrics.has_variant_qc:
            variant_df = qc_metrics.get_variant_metrics_df()
            variant_path = output_dir / f"{prefix}_variant_qc.csv"
            variant_df.to_csv(variant_path, index=False)
            saved_files['variant_qc'] = str(variant_path)
            logger.info(f"Saved variant QC metrics to {variant_path}")

        # Save MatrixTable if needed
        mt_path = output_dir / f"{prefix}_with_qc.mt"
        qc_metrics.mt.write(str(mt_path), overwrite=True)
        saved_files['matrix_table'] = str(mt_path)
        logger.info(f"Saved MatrixTable with QC annotations to {mt_path}")

        return saved_files

    except Exception as e:
        logger.error(f"Failed to save QC metrics: {e}")
        raise

