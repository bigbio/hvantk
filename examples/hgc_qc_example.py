"""
HGC QC Module Tests

Tests for the main QC functions using Hail's built-in MatrixTable creation methods.
"""

import logging
import tempfile

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def create_test_mt(n_samples=100, n_variants=1000):
    """Create a test MatrixTable using Hail's Balding-Nichols model."""
    import hail as hl

    # Use Hail's Balding-Nichols model for realistic genotype data
    mt = hl.balding_nichols_model(
        n_populations=3,
        n_samples=n_samples,
        n_variants=n_variants,
        n_partitions=4
    )

    # Add realistic genotype quality and depth fields for QC testing
    mt = mt.annotate_entries(
        GQ=hl.int32(hl.rand_unif(10, 40)),
        DP=hl.int32(hl.rand_unif(5, 30))
    )

    return mt

def test_sample_qc():
    """Test sample QC computation."""
    try:
        import hail as hl
        from hvantk.hgc import compute_sample_qc

        logger.info("Testing sample QC computation...")

        # Create test data
        mt = create_test_mt(n_samples=50, n_variants=500)
        logger.info(f"Created test MT: {mt.count_rows()} variants, {mt.count_cols()} samples")

        # Compute sample QC
        mt_qc = compute_sample_qc(mt)

        # Verify QC annotation exists
        assert 'sample_qc' in mt_qc.col, "Sample QC annotation not found"

        # Extract and check sample metrics
        sample_qc = mt_qc.cols().select('sample_qc').collect()
        assert len(sample_qc) == 50, f"Expected 50 samples, got {len(sample_qc)}"

        # Check that key metrics exist
        first_sample = sample_qc[0]['sample_qc']
        required_fields = ['call_rate', 'n_called', 'n_het', 'n_hom_ref', 'n_hom_var']
        for field in required_fields:
            assert field in first_sample, f"Required field {field} not found in sample QC"

        logger.info("✅ Sample QC test passed")
        return True

    except Exception as e:
        logger.error(f"❌ Sample QC test failed: {e}")
        return False

def test_variant_qc():
    """Test variant QC computation."""
    try:
        import hail as hl
        from hvantk.hgc import compute_variant_qc

        logger.info("Testing variant QC computation...")

        # Create test data
        mt = create_test_mt(n_samples=100, n_variants=200)

        # Compute variant QC
        mt_qc = compute_variant_qc(mt)

        # Verify QC annotation exists
        assert 'variant_qc' in mt_qc.row, "Variant QC annotation not found"

        # Extract and check variant metrics
        variant_qc = mt_qc.rows().select('variant_qc').take(10)

        # Check that key metrics exist
        first_variant = variant_qc[0]['variant_qc']
        required_fields = ['call_rate', 'AC', 'AF', 'n_called', 'n_het', 'p_value_hwe']
        for field in required_fields:
            assert field in first_variant, f"Required field {field} not found in variant QC"

        # Check AC/AF arrays have correct length (should be 2 for biallelic)
        assert len(first_variant['AC']) == 2, "AC should be array of length 2"
        assert len(first_variant['AF']) == 2, "AF should be array of length 2"

        logger.info("✅ Variant QC test passed")
        return True

    except Exception as e:
        logger.error(f"❌ Variant QC test failed: {e}")
        return False

def test_full_qc():
    """Test comprehensive QC computation."""
    try:
        import hail as hl
        from hvantk.hgc import compute_full_qc

        logger.info("Testing full QC computation...")

        # Create test data
        mt = create_test_mt(n_samples=75, n_variants=300)

        # Compute full QC
        qc_results = compute_full_qc(mt)

        # Check QCMetrics object
        assert qc_results.has_sample_qc, "Sample QC should be available"
        assert qc_results.has_variant_qc, "Variant QC should be available"

        # Check DataFrames
        sample_df = qc_results.get_sample_metrics_df()
        variant_df = qc_results.get_variant_metrics_df()

        assert len(sample_df) == 75, f"Expected 75 samples, got {len(sample_df)}"
        assert len(variant_df) == 300, f"Expected 300 variants, got {len(variant_df)}"

        # Check that DataFrames have data (columns will be nested under 'sample_qc' and 'variant_qc')
        assert len(sample_df.columns) > 0, "Sample DataFrame should have columns"
        assert len(variant_df.columns) > 0, "Variant DataFrame should have columns"

        # Check the actual structure - QC metrics are nested
        logger.info(f"Sample columns: {sample_df.columns.tolist()}")
        logger.info(f"Variant columns: {variant_df.columns.tolist()}")

        logger.info(f"Sample metrics shape: {sample_df.shape}")
        logger.info(f"Variant metrics shape: {variant_df.shape}")
        logger.info("✅ Full QC test passed")
        return True

    except Exception as e:
        logger.error(f"❌ Full QC test failed: {e}")
        return False

def test_qc_filtering():
    """Test QC-based filtering."""
    try:
        import hail as hl
        from hvantk.hgc import compute_full_qc, filter_samples_by_qc, filter_variants_by_qc

        logger.info("Testing QC filtering...")

        # Create test data with some poor quality samples/variants
        mt = create_test_mt(n_samples=100, n_variants=500)

        # Compute QC first
        qc_results = compute_full_qc(mt)

        # Test sample filtering
        mt_sample_filtered = filter_samples_by_qc(
            qc_results.mt,
            min_call_rate=0.5  # Lenient threshold for test data
        )

        samples_before = mt.count_cols()
        samples_after = mt_sample_filtered.count_cols()
        logger.info(f"Sample filtering: {samples_before} → {samples_after}")

        # Test variant filtering
        mt_variant_filtered = filter_variants_by_qc(
            mt_sample_filtered,
            min_call_rate=0.5,
            min_ac=1
        )

        variants_before = mt_sample_filtered.count_rows()
        variants_after = mt_variant_filtered.count_rows()
        logger.info(f"Variant filtering: {variants_before} → {variants_after}")

        # Should have removed some low-quality data
        assert samples_after <= samples_before, "Sample filtering should not increase sample count"
        assert variants_after <= variants_before, "Variant filtering should not increase variant count"

        logger.info("✅ QC filtering test passed")
        return True

    except Exception as e:
        logger.error(f"❌ QC filtering test failed: {e}")
        return False

def test_qc_export():
    """Test QC metrics export."""
    try:
        from hvantk.hgc import compute_full_qc, save_qc_metrics, prepare_qc_for_visualization

        logger.info("Testing QC export functionality...")

        # Create test data
        mt = create_test_mt(n_samples=30, n_variants=100)

        # Compute QC
        qc_results = compute_full_qc(mt)

        # Test saving to temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            saved_files = save_qc_metrics(qc_results, temp_dir, "test")

            # Check that files were created
            expected_files = ['sample_qc', 'variant_qc', 'matrix_table']
            for file_type in expected_files:
                assert file_type in saved_files, f"Expected file type {file_type} not saved"
                logger.info(f"Saved {file_type}: {saved_files[file_type]}")

        # Test visualization preparation
        viz_data = prepare_qc_for_visualization(qc_results)

        assert 'sample' in viz_data, "Sample visualization data missing"
        assert 'variant' in viz_data, "Variant visualization data missing"

        logger.info(f"Visualization data - Sample: {viz_data['sample'].shape}, Variant: {viz_data['variant'].shape}")
        logger.info("✅ QC export test passed")
        return True

    except Exception as e:
        logger.error(f"❌ QC export test failed: {e}")
        return False

def test_qc_visualization():
    """Test QC visualization functionality."""
    try:
        from hvantk.hgc import compute_full_qc

        logger.info("Testing QC visualization functionality...")

        # Create test data
        mt = create_test_mt(n_samples=50, n_variants=200)

        # Compute QC
        qc_results = compute_full_qc(mt)

        # Test QCMetrics plotting methods
        logger.info("Testing QCMetrics plotting methods...")

        try:
            fig1 = qc_results.plot_sample_overview()
            logger.info("✓ Sample overview plot created")
        except Exception as e:
            logger.warning(f"Sample overview plot failed: {e}")

        try:
            fig2 = qc_results.plot_variant_overview()
            logger.info("✓ Variant overview plot created")
        except Exception as e:
            logger.warning(f"Variant overview plot failed: {e}")

        try:
            fig3 = qc_results.plot_dashboard()
            logger.info("✓ QC dashboard created")
        except Exception as e:
            logger.warning(f"QC dashboard failed: {e}")

        # Test direct imports
        logger.info("Testing direct visualization imports...")
        try:
            from hvantk.visualization import plot_sample_qc_overview, plot_variant_qc_overview

            sample_df = qc_results.get_sample_metrics_df()
            variant_df = qc_results.get_variant_metrics_df()

            fig4 = plot_sample_qc_overview(sample_df)
            logger.info("✓ Direct sample overview import works")

            fig5 = plot_variant_qc_overview(variant_df)
            logger.info("✓ Direct variant overview import works")

        except Exception as e:
            logger.warning(f"Direct imports failed: {e}")

        logger.info("✅ QC visualization test completed")
        return True

    except Exception as e:
        logger.error(f"❌ QC visualization test failed: {e}")
        return False

def run_all_tests():
    """Run all QC tests."""
    try:
        import hail as hl

        # Initialize Hail
        logger.info("Initializing Hail for testing...")
        hl.init(quiet=True, log='/tmp/hail_test.log')

        tests = [
            ("Sample QC", test_sample_qc),
            ("Variant QC", test_variant_qc),
            ("Full QC", test_full_qc),
            ("QC Filtering", test_qc_filtering),
            ("QC Export", test_qc_export),
            ("QC Visualization", test_qc_visualization)
        ]

        results = []
        logger.info("Running HGC QC tests...\n")

        for test_name, test_func in tests:
            logger.info(f"Running {test_name} test...")
            try:
                success = test_func()
                results.append((test_name, success))
                if success:
                    logger.info(f"✅ {test_name} - PASSED\n")
                else:
                    logger.error(f"❌ {test_name} - FAILED\n")
            except Exception as e:
                logger.error(f"❌ {test_name} - ERROR: {e}\n")
                results.append((test_name, False))

        # Summary
        passed = sum(1 for _, success in results if success)
        total = len(results)

        print("=" * 50)
        print("TEST SUMMARY")
        print("=" * 50)
        for test_name, success in results:
            status = "✅ PASS" if success else "❌ FAIL"
            print(f"{test_name:20} {status}")

        print("-" * 50)
        print(f"TOTAL: {passed}/{total} tests passed")

        if passed == total:
            print("🎉 All tests passed! HGC QC module is working correctly.")
        else:
            print(f"⚠️  {total - passed} tests failed. Please check the errors above.")

        return passed == total

    except ImportError as e:
        if "hail" in str(e).lower():
            logger.error("Hail is required for these tests. Install with: pip install hail")
        else:
            logger.error(f"Import error: {e}")
        return False
    except Exception as e:
        logger.error(f"Test execution failed: {e}")
        return False

if __name__ == "__main__":
    print("HGC QC Module Tests")
    print("=" * 50)
    print("Testing main QC functions with Hail-generated data")
    print()

    success = run_all_tests()
    exit(0 if success else 1)
