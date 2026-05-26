"""PCA computation for ancestry inference.

This module provides functions to compute HWE-normalized PCA on merged
MatrixTables for ancestry inference, and to project new samples onto
existing PC space.
"""

import logging
from dataclasses import dataclass, field
from typing import List, Optional

import hail as hl
import numpy as np
import pandas as pd

from hvantk.algorithms.ancestry.constants import DEFAULT_N_PCS

logger = logging.getLogger(__name__)


@dataclass
class PCAResult:
    """Container for PCA results.

    Attributes
    ----------
    eigenvalues : List[float]
        Eigenvalues from PCA decomposition.
    scores : hl.Table
        Hail Table with PC scores per sample.
    loadings : hl.Table
        Hail Table with PC loadings per variant.
    n_variants : int
        Number of variants used in PCA.
    n_samples : int
        Number of samples in PCA.
    """

    eigenvalues: List[float]
    scores: hl.Table
    loadings: hl.Table
    n_variants: int
    n_samples: int
    _scores_df: Optional[pd.DataFrame] = field(default=None, repr=False)

    def get_scores_df(self) -> pd.DataFrame:
        """Convert scores to pandas DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame with sample IDs and PC scores (PC1, PC2, ..., PCn).
        """
        if self._scores_df is None:
            self._scores_df = self.scores.to_pandas()
        return self._scores_df

    def variance_explained(self) -> List[float]:
        """Return proportion of variance explained by each PC.

        Returns
        -------
        List[float]
            Proportion of variance explained by each principal component.
        """
        total = sum(self.eigenvalues)
        if total == 0:
            return [0.0] * len(self.eigenvalues)
        return [ev / total for ev in self.eigenvalues]

    def cumulative_variance(self) -> List[float]:
        """Return cumulative variance explained.

        Returns
        -------
        List[float]
            Cumulative proportion of variance explained up to each PC.
        """
        var_exp = self.variance_explained()
        return list(np.cumsum(var_exp))

    def get_n_pcs(self) -> int:
        """Return the number of PCs computed.

        Returns
        -------
        int
            Number of principal components.
        """
        return len(self.eigenvalues)


def compute_pca(
    mt: hl.MatrixTable,
    n_pcs: int = DEFAULT_N_PCS,
    compute_loadings: bool = True,
) -> PCAResult:
    """Compute HWE-normalized PCA on a MatrixTable.

    Performs principal component analysis using Hail's hwe_normalized_pca
    function, which accounts for allele frequency differences between
    populations (Patterson et al., 2006).

    Parameters
    ----------
    mt : hl.MatrixTable
        Input MatrixTable with GT field and variant_qc annotation.
        Should be filtered and LD-pruned before PCA.
    n_pcs : int, optional
        Number of principal components to compute. Default: 20.
    compute_loadings : bool, optional
        Whether to compute variant loadings (required for projection).
        Default: True.

    Returns
    -------
    PCAResult
        PCA results including eigenvalues, scores, and loadings.

    Raises
    ------
    ValueError
        If n_pcs <= 0 or if GT field is missing from MatrixTable.

    Notes
    -----
    The HWE-normalized PCA is equivalent to scaling genotypes by
    sqrt(p*(1-p)) where p is the allele frequency. This is standard
    for population genetics analysis.

    Example
    -------
    >>> from hvantk.algorithms.ancestry.pca import compute_pca
    >>> pca_result = compute_pca(filtered_mt, n_pcs=20)
    >>> print(f"Variance explained by PC1: {pca_result.variance_explained()[0]:.2%}")
    """
    # Validate input
    if n_pcs <= 0:
        raise ValueError(f"n_pcs must be positive, got {n_pcs}")

    if "GT" not in mt.entry:
        raise ValueError("GT field not found in MatrixTable entry fields")

    # Get counts for result object
    n_variants = mt.count_rows()
    n_samples = mt.count_cols()

    logger.info(
        f"Computing PCA with {n_pcs} PCs on {n_variants} variants, {n_samples} samples"
    )

    # Check for sufficient data to compute requested number of PCs
    # Mathematically, we can compute at most min(n_variants, n_samples) - 1 PCs
    min_dim = min(n_variants, n_samples)
    if min_dim <= n_pcs:
        raise ValueError(
            f"Insufficient data for PCA: need at least {n_pcs + 1} variants AND "
            f"{n_pcs + 1} samples to compute {n_pcs} PCs, but got {n_variants} "
            f"variants and {n_samples} samples. Maximum PCs possible: {min_dim - 1}"
        )

    # Warn if variant count is below recommended threshold for ancestry inference
    # gnomAD uses ~200k LD-pruned variants; research shows minimum ~100-1000 for
    # continental-level ancestry differentiation
    if n_variants < 1000:
        logger.warning(
            f"Only {n_variants} variants available for PCA. Ancestry inference may be "
            f"unreliable with fewer than 1,000 variants. gnomAD uses ~200k LD-pruned "
            f"variants for robust population structure analysis."
        )
    elif n_variants < 10000:
        logger.warning(
            f"Only {n_variants} variants available for PCA. For robust ancestry "
            f"inference, consider using 10,000+ LD-pruned variants."
        )

    # Warn if sample size is very small
    if n_samples < 50:
        logger.warning(
            f"Only {n_samples} samples available for PCA. Small sample sizes may "
            f"affect robustness of population structure analysis."
        )

    # Ensure n_pcs doesn't exceed minimum dimension
    max_pcs = min(n_variants, n_samples) - 1
    if n_pcs > max_pcs:
        logger.warning(
            f"Requested {n_pcs} PCs but maximum is {max_pcs}. Reducing to {max_pcs}."
        )
        n_pcs = max_pcs

    # Compute HWE-normalized PCA
    # Note: hwe_normalized_pca expects a call expression (GT), not n_alt_alleles()
    eigenvalues, scores_ht, loadings_ht = hl.hwe_normalized_pca(
        mt.GT,
        k=n_pcs,
        compute_loadings=compute_loadings,
    )

    # Annotate scores with individual PC columns for convenience
    pc_annotations = {f"PC{i + 1}": scores_ht.scores[i] for i in range(n_pcs)}
    scores_ht = scores_ht.annotate(**pc_annotations)

    # If loadings were computed, annotate with AF for future projection
    if compute_loadings and "variant_qc" in mt.row:
        loadings_ht = loadings_ht.annotate(
            af=mt.rows()[loadings_ht.key].variant_qc.AF[1]
        )

    # Log variance explained
    total_var = sum(eigenvalues)
    if total_var > 0:
        var_explained = [ev / total_var for ev in eigenvalues[:5]]
        logger.info(
            f"Variance explained by first 5 PCs: "
            f"{', '.join(f'{v:.2%}' for v in var_explained)}"
        )
        logger.info(f"Cumulative variance (PC1-5): {sum(var_explained):.2%}")

    return PCAResult(
        eigenvalues=list(eigenvalues),
        scores=scores_ht,
        loadings=loadings_ht if compute_loadings else None,
        n_variants=n_variants,
        n_samples=n_samples,
    )


def project_samples(
    mt: hl.MatrixTable,
    loadings: hl.Table,
    n_pcs: int = DEFAULT_N_PCS,
) -> hl.Table:
    """Project new samples onto existing PC space using loadings.

    Uses the variant loadings from a previous PCA to project new samples
    without recomputing the full PCA. This is useful for adding new samples
    to an existing ancestry model.

    Parameters
    ----------
    mt : hl.MatrixTable
        MatrixTable with new samples to project.
    loadings : hl.Table
        Variant loadings from previous PCA (from PCAResult.loadings).
        Must include the 'af' field for mean-centering.
    n_pcs : int, optional
        Number of PCs to project. Must not exceed loadings dimensions.
        Default: 20.

    Returns
    -------
    hl.Table
        Sample scores in projected PC space, with columns PC1, PC2, ..., PCn.

    Raises
    ------
    ValueError
        If loadings table is missing required fields or n_pcs is invalid.

    Notes
    -----
    The projection formula is:
        score[k] = sum_i ((g_i - 2*af_i) * loading_i[k])

    where g_i is the genotype (0, 1, 2), af_i is the allele frequency,
    and loading_i[k] is the k-th loading for variant i.

    Example
    -------
    >>> # Train PCA on reference
    >>> pca_result = compute_pca(reference_mt)
    >>> # Project new samples
    >>> new_scores = project_samples(new_mt, pca_result.loadings)
    """
    # Validate loadings
    if loadings is None:
        raise ValueError("loadings table is required for projection")

    if "loadings" not in loadings.row:
        raise ValueError("loadings field not found in loadings table")

    if "af" not in loadings.row:
        raise ValueError(
            "af field not found in loadings table. "
            "Ensure loadings were computed with variant_qc annotation."
        )

    if n_pcs <= 0:
        raise ValueError(f"n_pcs must be positive, got {n_pcs}")

    logger.info(f"Projecting {mt.count_cols()} samples onto {n_pcs} PCs")

    # Filter MT to variants in loadings
    mt = mt.semi_join_rows(loadings)
    n_overlap = mt.count_rows()
    logger.info(f"Projecting using {n_overlap} overlapping variants")

    if n_overlap == 0:
        raise ValueError("No overlapping variants between input MT and loadings")

    # Annotate MT with loadings
    mt = mt.annotate_rows(
        loadings=loadings[mt.row_key].loadings,
        af=loadings[mt.row_key].af,
    )

    # Compute projection: centered genotype * loadings
    # Use hl.pc_project which does this efficiently
    # Note: pc_project expects a call expression (GT), not n_alt_alleles()
    # Use the annotated row fields to avoid expression source mismatch
    scores_ht = hl.experimental.pc_project(
        mt.GT,
        mt.loadings,
        mt.af,
    )

    # Annotate with individual PC columns
    n_computed_pcs = min(n_pcs, len(scores_ht.scores.take(1)[0]))
    pc_annotations = {f"PC{i + 1}": scores_ht.scores[i] for i in range(n_computed_pcs)}
    scores_ht = scores_ht.annotate(**pc_annotations)

    return scores_ht
