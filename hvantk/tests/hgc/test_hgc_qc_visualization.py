"""
Test suite for HGC QC computation and HTML report generation.

This module contains tests to ensure QC computation and the HTML QC report
path work correctly. Tests create small MatrixTables, compute QC metrics,
and validate the resulting data and generated report.
"""

import pytest
import tempfile
from pathlib import Path

import hail as hl


@pytest.fixture(scope="session")
def test_mt():
    """Create a test MatrixTable once per session and reuse across all tests.

    This fixture creates a small test MatrixTable with realistic QC structure
    using Hail's Balding-Nichols model. Since MatrixTable creation is expensive,
    we use session scope to generate it once and share it across all tests.
    """
    # Create realistic test data using Balding-Nichols model
    mt = hl.balding_nichols_model(
        n_populations=2,  # Two populations for genetic diversity
        n_samples=25,  # Small but sufficient sample size
        n_variants=150,  # Enough variants for meaningful QC
        n_partitions=2,
    )

    # Add realistic entry fields that QC will analyze
    mt = mt.annotate_entries(
        GQ=hl.int32(hl.rand_unif(10, 40)),  # Genotype quality 10-40
        DP=hl.int32(hl.rand_unif(5, 35)),  # Depth 5-35
        AD=hl.array(
            [
                hl.int32(hl.rand_unif(0, 15)),  # Reference allele depth
                hl.int32(hl.rand_unif(0, 15)),  # Alternate allele depth
            ]
        ),
    )

    # Add some sample metadata for more realistic QC
    mt = mt.annotate_cols(population=hl.if_else(mt.sample_idx < 12, "POP1", "POP2"))

    return mt


def test_create_matrixtable_and_compute_qc(test_mt):
    """Test 1: Create MatrixTable and compute comprehensive QC metrics."""
    from hvantk.algorithms.hgc import compute_full_qc

    # Use the pre-created test MatrixTable from fixture
    mt = test_mt

    # Verify MatrixTable structure
    assert mt.count_cols() == 25, "Should have 25 samples"
    assert mt.count_rows() == 150, "Should have 150 variants"

    # Compute comprehensive QC
    qc_results = compute_full_qc(mt)

    # Verify QC computation
    assert qc_results.has_sample_qc, "Should have sample QC"
    assert qc_results.has_variant_qc, "Should have variant QC"

    # Check sample QC metrics
    sample_df = qc_results.get_sample_metrics_df()
    assert len(sample_df) == 25, "Sample QC should have 25 rows"

    # Verify key sample QC columns exist
    expected_sample_cols = [
        "sample_qc.call_rate",
        "sample_qc.n_called",
        "sample_qc.n_het",
        "sample_qc.n_hom_var",
        "sample_qc.r_ti_tv",
        "sample_qc.dp_stats.mean",
    ]
    for col in expected_sample_cols:
        assert col in sample_df.columns, f"Missing sample QC column: {col}"

    # Check variant QC metrics
    variant_df = qc_results.get_variant_metrics_df()
    assert len(variant_df) == 150, "Variant QC should have 150 rows"

    # Verify key variant QC columns exist
    expected_variant_cols = [
        "variant_qc.call_rate",
        "variant_qc.AF",
        "variant_qc.AC",
        "variant_qc.n_het",
        "variant_qc.p_value_hwe",
    ]
    for col in expected_variant_cols:
        assert col in variant_df.columns, f"Missing variant QC column: {col}"

    # Verify QC values are reasonable
    call_rates = sample_df["sample_qc.call_rate"]
    assert call_rates.min() >= 0.0, "Call rates should be >= 0"
    assert call_rates.max() <= 1.0, "Call rates should be <= 1"

    return qc_results


def test_html_qc_report_generation(test_mt):
    """Test 2: Generate comprehensive HTML QC reports."""
    from hvantk.algorithms.hgc import compute_full_qc
    from hvantk.algorithms.visualization.qc_report import generate_qc_report

    # Use the pre-created test MatrixTable from fixture
    mt = test_mt
    qc_results = compute_full_qc(mt)

    report_title = "Test QC Report"

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Test 1: Generate HTML report with all plots
        report_path = Path(tmp_dir) / "qc_report.html"
        result_path = generate_qc_report(qc_results, report_path, title=report_title)

        # Verify report was created
        assert result_path.exists(), "HTML report should be created"
        assert (
            result_path.stat().st_size > 50000
        ), "Report should be substantial (>50KB)"

        # Check report contains expected content
        report_content = result_path.read_text()
        assert (
            report_title in report_content
        ), "Report should contain the specified title"

        # Test 2: Generate report through QCMetrics
        report_path2 = Path(tmp_dir) / "qcmetrics_report.html"
        result_path2 = qc_results.generate_html_report(
            report_path2, title="QCMetrics Generated Report"
        )

        assert result_path2.exists(), "QCMetrics HTML report should be created"
        assert (
            result_path2.stat().st_size > 50000
        ), "QCMetrics report should be substantial"

        # Verify the title parameter works through QCMetrics wrapper
        report_content2 = result_path2.read_text()
        assert (
            "QCMetrics Generated Report" in report_content2
        ), "Report should contain the QCMetrics title"


def test_qc_data_validation(test_mt):
    """Test 3: Validate QC data quality and structure."""
    from hvantk.algorithms.hgc import compute_full_qc

    # Use the pre-created test MatrixTable from fixture
    mt = test_mt
    qc_results = compute_full_qc(mt)

    # Test sample QC data quality
    sample_df = qc_results.get_sample_metrics_df()

    # Check call rates are in valid range
    call_rates = sample_df["sample_qc.call_rate"]
    assert (call_rates >= 0).all(), "All sample call rates should be >= 0"
    assert (call_rates <= 1).all(), "All sample call rates should be <= 1"

    # Test variant QC data quality
    variant_df = qc_results.get_variant_metrics_df()

    # Check variant call rates
    var_call_rates = variant_df["variant_qc.call_rate"]
    assert (var_call_rates >= 0).all(), "All variant call rates should be >= 0"
    assert (var_call_rates <= 1).all(), "All variant call rates should be <= 1"

    # Check allele frequencies
    if "variant_qc.AF" in variant_df.columns:
        # Handle AF as array or float
        af_data = variant_df["variant_qc.AF"]
        if hasattr(af_data.iloc[0], "__len__") and len(af_data.iloc[0]) > 1:
            # AF is array format, check alternate allele frequency
            alt_afs = [af[1] if len(af) > 1 else 0 for af in af_data]
            assert all(
                0 <= af <= 1 for af in alt_afs
            ), "Allele frequencies should be in [0,1]"
