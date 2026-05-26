"""Variant filtering for ancestry inference.

This module provides functions for filtering variants to a high-quality,
informative subset suitable for ancestry analysis. The filtering pipeline
includes:

1. Autosomal chromosomes only
2. Biallelic SNPs
3. Allele frequency and call rate filters
4. Optional Hardy-Weinberg equilibrium filter
5. LD pruning for independent variants

Functions
---------
filter_to_autosomes
    Filter MatrixTable to autosomal chromosomes only.
filter_to_biallelic_snps
    Filter to biallelic SNP variants.
filter_variants_for_ancestry
    Apply ancestry-specific variant filters.
ld_prune_variants
    LD prune variants using Hail's ld_prune.
"""

import logging
from typing import Optional

import hail as hl

from hvantk.algorithms.ancestry.constants import (
    DEFAULT_MIN_AF,
    DEFAULT_MAX_AF,
    DEFAULT_MIN_CALL_RATE,
    DEFAULT_HWE_P,
    DEFAULT_LD_R2,
    DEFAULT_LD_WINDOW,
)

logger = logging.getLogger(__name__)


def filter_to_autosomes(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter MatrixTable to autosomal chromosomes only.

    Removes variants on sex chromosomes (X, Y) and mitochondrial DNA to avoid
    complications from different ploidy and inheritance patterns.

    Parameters
    ----------
    mt : hl.MatrixTable
        Input MatrixTable.

    Returns
    -------
    hl.MatrixTable
        MatrixTable filtered to autosomal variants only.

    Examples
    --------
    >>> mt_auto = filter_to_autosomes(mt)
    >>> # Only chr1-22 remain
    """
    n_before = mt.count_rows()

    mt = mt.filter_rows(mt.locus.in_autosome())

    n_after = mt.count_rows()
    logger.info(
        f"Autosome filter: {n_before:,} -> {n_after:,} variants "
        f"({n_before - n_after:,} removed)"
    )

    return mt


def filter_to_biallelic_snps(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter to biallelic SNP variants.

    Removes multi-allelic variants and indels, keeping only biallelic
    single nucleotide polymorphisms which are more reliable for ancestry
    analysis.

    Parameters
    ----------
    mt : hl.MatrixTable
        Input MatrixTable.

    Returns
    -------
    hl.MatrixTable
        MatrixTable filtered to biallelic SNPs only.

    Examples
    --------
    >>> mt_snps = filter_to_biallelic_snps(mt)
    """
    n_before = mt.count_rows()

    # Filter to exactly 2 alleles (ref + 1 alt)
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)

    # Filter to SNPs (both ref and alt are single nucleotides)
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))

    n_after = mt.count_rows()
    logger.info(
        f"Biallelic SNP filter: {n_before:,} -> {n_after:,} variants "
        f"({n_before - n_after:,} removed)"
    )

    return mt


def compute_variant_qc(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Compute variant QC metrics using Hail's built-in function.

    Adds `variant_qc` row annotation with metrics including call_rate, AF,
    n_called, n_not_called, etc.

    Parameters
    ----------
    mt : hl.MatrixTable
        Input MatrixTable with GT entry field.

    Returns
    -------
    hl.MatrixTable
        MatrixTable with variant_qc row annotation.

    Notes
    -----
    If variant_qc already exists, this function returns the MT unchanged
    to avoid redundant computation.
    """
    if "variant_qc" in mt.row:
        logger.info("variant_qc already computed, skipping")
        return mt

    logger.info("Computing variant QC metrics")
    mt = hl.variant_qc(mt)

    return mt


def filter_variants_for_ancestry(
    mt: hl.MatrixTable,
    min_af: float = DEFAULT_MIN_AF,
    max_af: float = DEFAULT_MAX_AF,
    min_call_rate: float = DEFAULT_MIN_CALL_RATE,
    apply_hwe_filter: bool = False,
    hwe_p_threshold: float = DEFAULT_HWE_P,
) -> hl.MatrixTable:
    """
    Apply ancestry-specific variant filters.

    Filters variants to a high-quality, informative subset suitable for
    ancestry analysis. The filtering steps are:
    1. Filter to autosomes
    2. Filter to biallelic SNPs
    3. Compute variant QC metrics
    4. Apply call rate filter
    5. Apply allele frequency filter
    6. Optionally apply HWE filter

    Parameters
    ----------
    mt : hl.MatrixTable
        Input MatrixTable (typically merged query + reference).
    min_af : float, optional
        Minimum allele frequency. Default is 0.01.
    max_af : float, optional
        Maximum allele frequency. Default is 0.99.
    min_call_rate : float, optional
        Minimum call rate. Default is 0.98.
    apply_hwe_filter : bool, optional
        Whether to apply Hardy-Weinberg equilibrium filter. Default is False.
    hwe_p_threshold : float, optional
        HWE p-value threshold; variants with p < threshold are removed.
        Default is 1e-6.

    Returns
    -------
    hl.MatrixTable
        Filtered MatrixTable with variant_qc annotation.

    Raises
    ------
    ValueError
        If min_af >= max_af or if parameters are out of valid range.

    Examples
    --------
    >>> mt_filtered = filter_variants_for_ancestry(
    ...     mt,
    ...     min_af=0.05,
    ...     min_call_rate=0.99,
    ... )
    """
    # Validate parameters
    if min_af < 0 or min_af >= 1:
        raise ValueError(f"min_af must be in [0, 1), got {min_af}")
    if max_af <= 0 or max_af > 1:
        raise ValueError(f"max_af must be in (0, 1], got {max_af}")
    if min_af >= max_af:
        raise ValueError(f"min_af ({min_af}) must be less than max_af ({max_af})")
    if min_call_rate < 0 or min_call_rate > 1:
        raise ValueError(f"min_call_rate must be in [0, 1], got {min_call_rate}")

    n_initial = mt.count_rows()
    logger.info(f"Starting variant filtering with {n_initial:,} variants")

    # Step 1: Filter to autosomes
    mt = filter_to_autosomes(mt)

    # Step 2: Filter to biallelic SNPs
    mt = filter_to_biallelic_snps(mt)

    # Step 3: Compute variant QC
    mt = compute_variant_qc(mt)

    # Step 4: Apply call rate filter
    n_before = mt.count_rows()
    mt = mt.filter_rows(mt.variant_qc.call_rate >= min_call_rate)
    n_after = mt.count_rows()
    logger.info(
        f"Call rate filter (>= {min_call_rate}): {n_before:,} -> {n_after:,} variants"
    )

    # Step 5: Apply AF filter (on alt allele frequency)
    n_before = n_after
    # AF is an array; AF[1] is the alt allele frequency
    mt = mt.filter_rows(
        (mt.variant_qc.AF[1] >= min_af) & (mt.variant_qc.AF[1] <= max_af)
    )
    n_after = mt.count_rows()
    logger.info(
        f"AF filter ({min_af} <= AF <= {max_af}): {n_before:,} -> {n_after:,} variants"
    )

    # Step 6: Optional HWE filter
    if apply_hwe_filter:
        n_before = n_after
        # HWE p-value is in variant_qc.p_value_hwe
        mt = mt.filter_rows(mt.variant_qc.p_value_hwe >= hwe_p_threshold)
        n_after = mt.count_rows()
        logger.info(
            f"HWE filter (p >= {hwe_p_threshold}): {n_before:,} -> {n_after:,} variants"
        )

    logger.info(
        f"Variant filtering complete: {n_initial:,} -> {n_after:,} variants "
        f"({100 * n_after / n_initial:.1f}% retained)"
    )

    return mt


def ld_prune_variants(
    mt: hl.MatrixTable,
    r2: float = DEFAULT_LD_R2,
    bp_window_size: int = DEFAULT_LD_WINDOW,
    memory_per_core: int = 512,
    block_size: Optional[int] = None,
) -> hl.MatrixTable:
    """
    LD prune variants using Hail's ld_prune.

    Removes correlated variants to ensure independent signals for PCA.
    Uses a sliding window approach to compute pairwise LD and prunes
    variants exceeding the r² threshold.

    Parameters
    ----------
    mt : hl.MatrixTable
        Input MatrixTable with GT field.
    r2 : float, optional
        LD r² threshold; variant pairs with r² >= threshold have one removed.
        Default is 0.2.
    bp_window_size : int, optional
        Window size in base pairs for LD computation. Default is 500,000 (500 kb).
    memory_per_core : int, optional
        Memory per core in MB for LD matrix computation. Default is 512.
    block_size : int, optional
        Block size for LD matrix computation. If None, uses Hail's default.

    Returns
    -------
    hl.MatrixTable
        MatrixTable filtered to LD-independent variants.

    Notes
    -----
    LD pruning can be computationally intensive for large datasets.
    Consider checkpointing before this step for large analyses.

    Examples
    --------
    >>> mt_pruned = ld_prune_variants(mt, r2=0.1, bp_window_size=250_000)
    """
    n_before = mt.count_rows()
    logger.info(
        f"Starting LD pruning: {n_before:,} variants, "
        f"r²={r2}, window={bp_window_size:,} bp"
    )

    # Prepare kwargs for ld_prune
    ld_kwargs = {
        "call_expr": mt.GT,
        "r2": r2,
        "bp_window_size": bp_window_size,
        "memory_per_core": memory_per_core,
    }

    if block_size is not None:
        ld_kwargs["block_size"] = block_size

    # Compute pruned variant table
    logger.info("Computing LD pruning (this may take a while for large datasets)")
    pruned_variant_table = hl.ld_prune(**ld_kwargs)

    # Filter MT to pruned variants
    mt = mt.filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))

    n_after = mt.count_rows()
    logger.info(
        f"LD pruning complete: {n_before:,} -> {n_after:,} variants "
        f"({n_before - n_after:,} removed, {100 * n_after / n_before:.1f}% retained)"
    )

    return mt


def prepare_ancestry_variants(
    mt: hl.MatrixTable,
    min_af: float = DEFAULT_MIN_AF,
    max_af: float = DEFAULT_MAX_AF,
    min_call_rate: float = DEFAULT_MIN_CALL_RATE,
    apply_hwe_filter: bool = False,
    hwe_p_threshold: float = DEFAULT_HWE_P,
    ld_r2: float = DEFAULT_LD_R2,
    ld_window: int = DEFAULT_LD_WINDOW,
    skip_ld_pruning: bool = False,
) -> hl.MatrixTable:
    """
    Prepare variants for ancestry analysis with all filtering steps.

    Convenience function that applies variant filtering and LD pruning
    in a single call.

    Parameters
    ----------
    mt : hl.MatrixTable
        Input MatrixTable.
    min_af : float, optional
        Minimum allele frequency. Default is 0.01.
    max_af : float, optional
        Maximum allele frequency. Default is 0.99.
    min_call_rate : float, optional
        Minimum call rate. Default is 0.98.
    apply_hwe_filter : bool, optional
        Whether to apply HWE filter. Default is False.
    hwe_p_threshold : float, optional
        HWE p-value threshold. Default is 1e-6.
    ld_r2 : float, optional
        LD r² threshold for pruning. Default is 0.2.
    ld_window : int, optional
        LD window size in bp. Default is 500,000.
    skip_ld_pruning : bool, optional
        If True, skip LD pruning step. Default is False.

    Returns
    -------
    hl.MatrixTable
        MatrixTable with filtered and optionally LD-pruned variants.

    Examples
    --------
    >>> mt_prepared = prepare_ancestry_variants(
    ...     mt,
    ...     min_af=0.05,
    ...     ld_r2=0.1,
    ... )
    """
    # Apply standard filters
    mt = filter_variants_for_ancestry(
        mt=mt,
        min_af=min_af,
        max_af=max_af,
        min_call_rate=min_call_rate,
        apply_hwe_filter=apply_hwe_filter,
        hwe_p_threshold=hwe_p_threshold,
    )

    # Apply LD pruning unless skipped
    if not skip_ld_pruning:
        mt = ld_prune_variants(
            mt=mt,
            r2=ld_r2,
            bp_window_size=ld_window,
        )

    return mt
