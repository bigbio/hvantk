"""Test fixtures for ancestry inference module tests.

Provides synthetic MatrixTables for testing ancestry functionality without
requiring real genetic data.
"""

import pytest


@pytest.fixture(scope="module")
def synthetic_query_mt(hail_session, tmp_path_factory):
    """
    Create synthetic query MatrixTable for testing.

    Generates 50 samples with 1000 variants using Hail's Balding-Nichols
    population genetics simulator. This represents a cohort with unknown
    ancestry to be classified.

    Returns
    -------
    hl.MatrixTable
        Synthetic query MatrixTable.
    """
    import hail as hl

    tmp_dir = tmp_path_factory.mktemp("ancestry_test")
    mt_path = str(tmp_dir / "query.mt")

    # Generate synthetic genotypes from a single population
    mt = hl.balding_nichols_model(
        n_populations=1,
        n_samples=50,
        n_variants=1000,
        n_partitions=4,
    )

    # Add sample IDs (convert int to string)
    mt = mt.annotate_cols(s=hl.str("query_") + hl.str(mt.sample_idx))
    mt = mt.key_cols_by("s")

    # Ensure we have standard locus/alleles keys
    mt = mt.key_rows_by("locus", "alleles")

    mt.write(mt_path, overwrite=True)
    return hl.read_matrix_table(mt_path)


@pytest.fixture(scope="module")
def synthetic_reference_mt(hail_session, tmp_path_factory):
    """
    Create synthetic reference MatrixTable with ancestry labels for testing.

    Generates 90 samples from 3 populations (30 each: EUR, AFR, EAS) with
    1000 variants using Hail's Balding-Nichols simulator with FST values
    that create distinct population clusters.

    Returns
    -------
    hl.MatrixTable
        Synthetic reference MatrixTable with ancestry annotations.
    """
    import hail as hl

    tmp_dir = tmp_path_factory.mktemp("ancestry_test")
    mt_path = str(tmp_dir / "reference.mt")

    # Generate 3 populations with distinct FST values
    mt = hl.balding_nichols_model(
        n_populations=3,
        n_samples=90,  # 30 per population
        n_variants=1000,
        pop_dist=[1 / 3, 1 / 3, 1 / 3],
        fst=[0.1, 0.1, 0.1],
        n_partitions=4,
    )

    # Create population labels based on sample index
    # pop_id is 0, 1, or 2 based on which population the sample belongs to
    pop_labels = hl.literal(["EUR", "AFR", "EAS"])

    mt = mt.annotate_cols(
        s=hl.str("ref_") + hl.str(mt.sample_idx),
        ancestry=pop_labels[mt.pop],
    )
    mt = mt.key_cols_by("s")

    # Ensure we have standard locus/alleles keys
    mt = mt.key_rows_by("locus", "alleles")

    mt.write(mt_path, overwrite=True)
    return hl.read_matrix_table(mt_path)


@pytest.fixture(scope="module")
def synthetic_reference_mt_single_pop(hail_session, tmp_path_factory):
    """
    Create synthetic reference MatrixTable with only one population.

    Used for testing error handling when reference has insufficient
    population diversity.

    Returns
    -------
    hl.MatrixTable
        Synthetic reference MatrixTable with single ancestry.
    """
    import hail as hl

    tmp_dir = tmp_path_factory.mktemp("ancestry_test")
    mt_path = str(tmp_dir / "reference_single.mt")

    mt = hl.balding_nichols_model(
        n_populations=1,
        n_samples=30,
        n_variants=1000,
        n_partitions=4,
    )

    mt = mt.annotate_cols(
        s=hl.str("ref_") + hl.str(mt.sample_idx),
        ancestry=hl.literal("EUR"),
    )
    mt = mt.key_cols_by("s")
    mt = mt.key_rows_by("locus", "alleles")

    mt.write(mt_path, overwrite=True)
    return hl.read_matrix_table(mt_path)


@pytest.fixture(scope="module")
def synthetic_incompatible_mt(hail_session, tmp_path_factory):
    """
    Create synthetic MatrixTable with different row key structure.

    Used for testing validation error handling.

    Returns
    -------
    hl.MatrixTable
        MatrixTable with non-standard key structure.
    """
    import hail as hl

    tmp_dir = tmp_path_factory.mktemp("ancestry_test")
    mt_path = str(tmp_dir / "incompatible.mt")

    mt = hl.balding_nichols_model(
        n_populations=1,
        n_samples=10,
        n_variants=100,
        n_partitions=2,
    )

    # Key by locus only (not alleles) to create incompatibility
    mt = mt.key_rows_by("locus")
    mt = mt.annotate_cols(s=hl.str(mt.sample_idx))
    mt = mt.key_cols_by("s")

    mt.write(mt_path, overwrite=True)
    return hl.read_matrix_table(mt_path)


@pytest.fixture(scope="module")
def small_merged_mt(hail_session, synthetic_query_mt, synthetic_reference_mt):
    """
    Pre-merged MatrixTable for filter testing.

    Merges query and reference MTs to provide a fixture for testing
    filtering functions without re-running the merge each time.

    Returns
    -------
    hl.MatrixTable
        Merged MatrixTable.
    """
    from hvantk.ancestry.merge import merge_matrixtables

    return merge_matrixtables(
        query_mt=synthetic_query_mt,
        reference_mt=synthetic_reference_mt,
        ancestry_col="ancestry",
        validate=True,
        min_shared_variants=100,  # Lower threshold for testing
    )
