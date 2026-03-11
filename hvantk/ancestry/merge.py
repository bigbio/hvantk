"""MatrixTable merging utilities for ancestry inference.

This module provides functions for validating and merging query and reference
MatrixTables for ancestry analysis. The merge operation performs an inner join
on variants, keeping only variants present in both datasets.

Functions
---------
validate_matrixtable_compatibility
    Check that two MatrixTables are compatible for merging.
get_shared_variants_stats
    Count variants unique to each MT and shared between them.
merge_matrixtables
    Merge query and reference MatrixTables by shared variants.
"""

import logging
from typing import Optional, Tuple

import hail as hl

from hvantk.ancestry.constants import (
    ANCESTRY_COL,
    SOURCE_COL,
    KNOWN_ANCESTRY_COL,
    MIN_SHARED_VARIANTS,
    WARN_SHARED_VARIANTS,
)

logger = logging.getLogger(__name__)


def validate_matrixtable_compatibility(
    query_mt: hl.MatrixTable,
    reference_mt: hl.MatrixTable,
    ancestry_col: str = ANCESTRY_COL,
) -> None:
    """
    Validate that two MatrixTables are compatible for merging.

    Performs the following checks:
    - Same row key structure (locus, alleles)
    - Same reference genome
    - GT field present in both MatrixTables
    - Ancestry column present in reference MatrixTable

    Parameters
    ----------
    query_mt : hl.MatrixTable
        Query cohort MatrixTable.
    reference_mt : hl.MatrixTable
        Reference panel MatrixTable.
    ancestry_col : str, optional
        Name of the column containing ancestry labels in reference_mt.
        Default is "ancestry".

    Raises
    ------
    ValueError
        If MatrixTables are incompatible for merging.

    Examples
    --------
    >>> validate_matrixtable_compatibility(query_mt, reference_mt)
    >>> # No exception means MTs are compatible
    """
    # Check row key structure
    query_row_key = list(query_mt.row_key.keys())
    ref_row_key = list(reference_mt.row_key.keys())

    if query_row_key != ref_row_key:
        raise ValueError(
            f"Row key mismatch: query has {query_row_key}, "
            f"reference has {ref_row_key}. "
            "Both must be keyed by (locus, alleles)."
        )

    # Check for expected key structure
    expected_keys = ["locus", "alleles"]
    if query_row_key != expected_keys:
        raise ValueError(
            f"Unexpected row key structure: {query_row_key}. "
            f"Expected {expected_keys}."
        )

    # Check reference genome compatibility
    query_locus_dtype = query_mt.locus.dtype
    ref_locus_dtype = reference_mt.locus.dtype

    if query_locus_dtype.reference_genome != ref_locus_dtype.reference_genome:
        raise ValueError(
            f"Reference genome mismatch: query uses "
            f"{query_locus_dtype.reference_genome.name}, "
            f"reference uses {ref_locus_dtype.reference_genome.name}."
        )

    # Check GT field presence
    query_entry_fields = list(query_mt.entry)
    ref_entry_fields = list(reference_mt.entry)

    if "GT" not in query_entry_fields and "LGT" not in query_entry_fields:
        raise ValueError(
            "Query MatrixTable must have GT or LGT entry field. "
            f"Found entry fields: {query_entry_fields}"
        )

    if "GT" not in ref_entry_fields and "LGT" not in ref_entry_fields:
        raise ValueError(
            "Reference MatrixTable must have GT or LGT entry field. "
            f"Found entry fields: {ref_entry_fields}"
        )

    # Check ancestry column presence in reference
    ref_col_fields = list(reference_mt.col)
    if ancestry_col not in ref_col_fields:
        raise ValueError(
            f"Ancestry column '{ancestry_col}' not found in reference MatrixTable. "
            f"Available column fields: {ref_col_fields}"
        )

    logger.info(
        f"MatrixTables validated: both use {query_locus_dtype.reference_genome.name}, "
        f"keyed by {query_row_key}"
    )


def get_shared_variants_stats(
    mt1: hl.MatrixTable,
    mt2: hl.MatrixTable,
) -> Tuple[int, int, int]:
    """
    Count variants unique to each MatrixTable and shared between them.

    Parameters
    ----------
    mt1 : hl.MatrixTable
        First MatrixTable.
    mt2 : hl.MatrixTable
        Second MatrixTable.

    Returns
    -------
    Tuple[int, int, int]
        Tuple of (mt1_only, mt2_only, shared) variant counts.

    Examples
    --------
    >>> mt1_only, mt2_only, shared = get_shared_variants_stats(query_mt, ref_mt)
    >>> print(f"Shared variants: {shared}")
    """
    # Get row keys as tables
    rows1 = mt1.rows().select()
    rows2 = mt2.rows().select()

    # Count totals
    n1 = rows1.count()
    n2 = rows2.count()

    # Count shared (inner join)
    shared_rows = rows1.semi_join(rows2)
    n_shared = shared_rows.count()

    # Calculate unique counts
    n1_only = n1 - n_shared
    n2_only = n2 - n_shared

    logger.info(
        f"Variant overlap: MT1={n1:,} ({n1_only:,} unique), "
        f"MT2={n2:,} ({n2_only:,} unique), shared={n_shared:,}"
    )

    return n1_only, n2_only, n_shared


def _normalize_gt_field(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Normalize genotype field to GT if LGT is present.

    Some MatrixTables (especially from VDS conversion) have LGT instead of GT.
    This function ensures a consistent GT field is present.

    Parameters
    ----------
    mt : hl.MatrixTable
        Input MatrixTable.

    Returns
    -------
    hl.MatrixTable
        MatrixTable with GT field.
    """
    if "GT" in mt.entry:
        return mt
    elif "LGT" in mt.entry:
        # Validate that LA field exists for lgt_to_gt conversion
        if "LA" not in mt.entry:
            raise ValueError(
                "MatrixTable has LGT entry field but is missing required LA (local alleles) field. "
                "LA is required for converting LGT to GT."
            )
        logger.info("Converting LGT to GT")
        return mt.annotate_entries(GT=hl.vds.lgt_to_gt(mt.LGT, mt.LA))
    else:
        raise ValueError("MatrixTable must have GT or LGT entry field")


def merge_matrixtables(
    query_mt: hl.MatrixTable,
    reference_mt: hl.MatrixTable,
    ancestry_col: str = ANCESTRY_COL,
    validate: bool = True,
    min_shared_variants: Optional[int] = None,
) -> hl.MatrixTable:
    """
    Merge query and reference MatrixTables by shared variants.

    Performs an inner join on variant keys, keeping only variants present in
    both MatrixTables. Samples from both datasets are combined via column
    union, with source and ancestry information annotated.

    Parameters
    ----------
    query_mt : hl.MatrixTable
        Query cohort with unknown ancestry.
    reference_mt : hl.MatrixTable
        Reference panel with known ancestry labels.
    ancestry_col : str, optional
        Column name containing ancestry labels in reference_mt.
        Default is "ancestry".
    validate : bool, optional
        Whether to validate compatibility before merging. Default is True.
    min_shared_variants : int, optional
        Minimum number of shared variants required. Default is MIN_SHARED_VARIANTS.
        Set to a lower value for testing purposes.

    Returns
    -------
    hl.MatrixTable
        Merged MatrixTable with additional column annotations:
        - _ancestry_source: 'query' or 'reference'
        - _known_ancestry: Original ancestry label (None for query samples)

    Raises
    ------
    ValueError
        If no shared variants found or MatrixTables are incompatible.

    Notes
    -----
    Uses inner join: only variants present in both MTs are retained.
    The merge preserves the genotype field (GT) and drops other entry fields
    to ensure compatibility.

    Examples
    --------
    >>> merged_mt = merge_matrixtables(
    ...     query_mt=hl.read_matrix_table("cohort.mt"),
    ...     reference_mt=hl.read_matrix_table("1kg.mt"),
    ...     ancestry_col="super_pop",
    ... )
    >>> print(f"Merged samples: {merged_mt.count_cols()}")
    """
    if validate:
        validate_matrixtable_compatibility(query_mt, reference_mt, ancestry_col)

    # Log input statistics
    n_query_samples = query_mt.count_cols()
    n_ref_samples = reference_mt.count_cols()
    logger.info(f"Query samples: {n_query_samples:,}")
    logger.info(f"Reference samples: {n_ref_samples:,}")

    # Get shared variant statistics
    query_only, ref_only, n_shared = get_shared_variants_stats(query_mt, reference_mt)

    # Set minimum threshold
    min_variants = (
        min_shared_variants if min_shared_variants is not None else MIN_SHARED_VARIANTS
    )

    # Validate sufficient shared variants
    if n_shared == 0:
        raise ValueError(
            "No shared variants between query and reference MatrixTables. "
            "Ensure both use the same reference genome and variant representation."
        )

    if n_shared < min_variants:
        raise ValueError(
            f"Only {n_shared:,} shared variants found; need at least "
            f"{min_variants:,} for reliable ancestry inference."
        )

    if n_shared < WARN_SHARED_VARIANTS:
        logger.warning(
            f"Only {n_shared:,} shared variants found. Results may be less reliable "
            f"with fewer than {WARN_SHARED_VARIANTS:,} variants."
        )

    # Get shared variant keys
    shared_variants = query_mt.rows().select().semi_join(reference_mt.rows().select())

    # Filter both MTs to shared variants
    logger.info(f"Filtering to {n_shared:,} shared variants")
    query_filtered = query_mt.semi_join_rows(shared_variants)
    ref_filtered = reference_mt.semi_join_rows(shared_variants)

    # Normalize GT field
    query_filtered = _normalize_gt_field(query_filtered)
    ref_filtered = _normalize_gt_field(ref_filtered)

    # Select only GT entry field for compatibility
    query_filtered = query_filtered.select_entries("GT")
    ref_filtered = ref_filtered.select_entries("GT")

    # Drop row annotations to avoid conflicts during union
    query_filtered = query_filtered.select_rows()
    ref_filtered = ref_filtered.select_rows()

    # Annotate samples with source and known ancestry
    query_annotated = query_filtered.annotate_cols(
        **{
            SOURCE_COL: "query",
            KNOWN_ANCESTRY_COL: hl.missing(hl.tstr),
        }
    )

    ref_annotated = ref_filtered.annotate_cols(
        **{
            SOURCE_COL: "reference",
            KNOWN_ANCESTRY_COL: ref_filtered[ancestry_col],
        }
    )

    # Drop other column annotations except the key and our new fields
    # First identify what columns to keep
    cols_to_keep = [SOURCE_COL, KNOWN_ANCESTRY_COL]

    # Select only our annotation columns (key is preserved automatically)
    query_annotated = query_annotated.select_cols(*cols_to_keep)
    ref_annotated = ref_annotated.select_cols(*cols_to_keep)

    # Union columns (combine samples)
    logger.info("Merging samples via column union")
    merged_mt = query_annotated.union_cols(ref_annotated)

    # Log final statistics
    n_merged_samples = merged_mt.count_cols()
    n_merged_variants = merged_mt.count_rows()
    logger.info(
        f"Merged MatrixTable: {n_merged_samples:,} samples, "
        f"{n_merged_variants:,} variants"
    )

    return merged_mt


def check_sample_overlap(
    query_mt: hl.MatrixTable,
    reference_mt: hl.MatrixTable,
) -> Tuple[int, hl.Table]:
    """
    Check for sample ID overlap between query and reference MatrixTables.

    Samples appearing in both datasets may indicate data leakage and should
    be handled appropriately (e.g., excluded from training).

    Parameters
    ----------
    query_mt : hl.MatrixTable
        Query cohort MatrixTable.
    reference_mt : hl.MatrixTable
        Reference panel MatrixTable.

    Returns
    -------
    Tuple[int, hl.Table]
        Tuple of (count, table) where count is the number of overlapping
        samples and table contains the overlapping sample IDs.

    Examples
    --------
    >>> n_overlap, overlap_samples = check_sample_overlap(query_mt, ref_mt)
    >>> if n_overlap > 0:
    ...     print(f"Warning: {n_overlap} samples appear in both datasets")
    """
    query_samples = query_mt.cols().select()
    ref_samples = reference_mt.cols().select()

    overlap = query_samples.semi_join(ref_samples)
    n_overlap = overlap.count()

    if n_overlap > 0:
        logger.warning(
            f"Found {n_overlap} sample(s) present in both query and reference. "
            "These may cause data leakage if used for training."
        )

    return n_overlap, overlap
