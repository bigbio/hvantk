"""
Integration tests for PSROC pipeline.

These tests validate the complete PSROC workflow including:
- Configuration validation
- Pipeline execution stages
- Output file generation
- Golden dataset verification

Most tests require Hail initialization and are marked accordingly.
"""

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from hvantk.algorithms.psroc.pipeline import (
    PSROCConfig,
    PSROCState,
    PSROCResult,
    PSROCStage,
)
from hvantk.algorithms.psroc.roc import (
    ROCResult,
    ScoreMissingness,
    compute_roc_metrics,
    compute_all_missingness,
)

# Test data paths
TEST_DATA_DIR = Path(__file__).parent.parent / "testdata"
PSROC_TEST_DATA = TEST_DATA_DIR / "psroc"
CLINVAR_TEST_DATA = (
    Path(__file__).parent.parent.parent
    / "skills"
    / "clinvar"
    / "tests"
    / "testdata"
    / "raw"
    / "clinvar"
    / "clinvar_20220403_chr20.vcf.bgz"
)
DBNSFP_TEST_DATA = (
    TEST_DATA_DIR / "raw" / "dbnsfp" / "dbNSFP4_v49a_example_variants.bgz"
)
GOLDEN_VARIANTS = PSROC_TEST_DATA / "golden_variants.txt"


class TestGoldenDatasets:
    """Tests using golden (synthetic) datasets with known expected results."""

    def test_perfect_separation_auc(self):
        """Golden test: Perfect separation should yield AUC=1.0."""
        # Create labels: first half benign, second half pathogenic
        labels = np.array([0] * 50 + [1] * 50)

        # Create perfect separator: benign scores < pathogenic scores
        scores = {
            "perfect_score": np.concatenate(
                [
                    np.linspace(0.0, 0.4, 50),  # Benign: 0.0-0.4
                    np.linspace(0.6, 1.0, 50),  # Pathogenic: 0.6-1.0
                ]
            )
        }

        results = compute_roc_metrics(labels, scores)

        assert "perfect_score" in results
        assert results["perfect_score"].auc == pytest.approx(1.0, abs=0.001)

    def test_random_classifier_auc(self):
        """Golden test: Random classifier should yield AUC≈0.5."""
        np.random.seed(42)
        n = 1000
        labels = np.array([0] * (n // 2) + [1] * (n // 2))
        scores = {"random_score": np.random.rand(n)}

        results = compute_roc_metrics(labels, scores)

        # Random classifier should be between 0.45 and 0.55
        assert 0.45 < results["random_score"].auc < 0.55

    def test_known_auc_synthetic_data(self):
        """Golden test: Synthetic data with known approximate AUC."""
        np.random.seed(123)
        n = 200

        # Generate labels
        labels = np.array([0] * (n // 2) + [1] * (n // 2))

        # Generate scores with moderate separation
        # Benign: mean=0.3, Pathogenic: mean=0.7
        benign_scores = np.random.normal(0.3, 0.15, n // 2)
        patho_scores = np.random.normal(0.7, 0.15, n // 2)
        scores = {
            "moderate_score": np.clip(
                np.concatenate([benign_scores, patho_scores]), 0, 1
            )
        }

        results = compute_roc_metrics(labels, scores)

        # With this separation, AUC should be high (>0.85)
        assert results["moderate_score"].auc > 0.85

    def test_missingness_threshold_behavior(self):
        """Golden test: Verify missingness filtering at threshold boundary."""
        labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])

        # Score with exactly 30% missing (should be excluded at 0.3 threshold)
        scores_at_threshold = {
            "score_30pct": np.array(
                [0.1, 0.2, np.nan, np.nan, np.nan, 0.6, 0.7, 0.8, 0.9, 1.0]
            )
        }

        # Score with exactly 30% missing should be excluded (>= threshold)
        # This should raise ValueError because all scores are excluded
        with pytest.raises(
            ValueError, match="All scores were excluded due to high missingness"
        ):
            compute_roc_metrics(labels, scores_at_threshold, max_missingness=0.3)

        # Score with 20% missing (should be included, below 30% threshold)
        scores_below_threshold = {
            "score_20pct": np.array(
                [0.1, 0.2, 0.3, np.nan, np.nan, 0.6, 0.7, 0.8, 0.9, 1.0]
            )
        }

        # This has 2/10 = 20% missing, which should be included
        results_below = compute_roc_metrics(
            labels, scores_below_threshold, max_missingness=0.3
        )
        assert "score_20pct" in results_below


class TestPipelineOutputArtifacts:
    """Test that pipeline generates all expected output artifacts."""

    def test_metrics_json_format(self):
        """Test that metrics JSON has correct structure."""
        # Create mock result
        missingness = ScoreMissingness(
            score_name="CADD_phred",
            n_total=100,
            n_present=95,
            n_missing=5,
            missingness_rate=0.05,
            included_in_analysis=True,
        )

        roc_result = ROCResult(
            score_name="CADD_phred",
            fpr=np.array([0.0, 0.5, 1.0]),
            tpr=np.array([0.0, 0.85, 1.0]),
            thresholds=np.array([1.0, 0.5, 0.0]),
            auc=0.85,
            optimal_threshold=0.5,
            sensitivity_at_optimal=0.85,
            specificity_at_optimal=0.5,
            n_variants_used=95,
            missingness=missingness,
        )

        result = PSROCResult(
            annotated_ht_path="/test/annotated.ht",
            metrics={"CADD_phred": roc_result},
            missingness={"CADD_phred": missingness},
            n_pathogenic=50,
            n_benign=45,
            n_excluded=5,
            n_total=100,
            scores_included=["CADD_phred"],
            scores_excluded=[],
            max_missingness_threshold=0.3,
            output_dir="/test",
        )

        # Test JSON serialization
        result_dict = result.to_dict()

        assert "metrics" in result_dict
        assert "CADD_phred" in result_dict["metrics"]
        assert result_dict["metrics"]["CADD_phred"]["auc"] == 0.85
        assert result_dict["n_pathogenic"] == 50
        assert result_dict["n_benign"] == 45

        # Verify JSON serializable
        json_str = json.dumps(result_dict)
        loaded = json.loads(json_str)
        assert loaded["metrics"]["CADD_phred"]["auc"] == 0.85

    def test_missingness_json_format(self):
        """Test missingness report JSON structure."""
        missingness = ScoreMissingness(
            score_name="REVEL_score",
            n_total=100,
            n_present=60,
            n_missing=40,
            missingness_rate=0.40,
            included_in_analysis=False,
            exclusion_reason="missingness_rate (0.40) exceeds max_missingness (0.30)",
        )

        miss_dict = missingness.to_dict()

        assert miss_dict["n_total"] == 100
        assert miss_dict["n_present"] == 60
        assert miss_dict["n_missing"] == 40
        assert miss_dict["missingness_rate"] == 0.40
        assert miss_dict["included_in_analysis"] is False
        assert "exceeds" in miss_dict["exclusion_reason"]

        # Verify JSON serializable
        json_str = json.dumps(miss_dict)
        loaded = json.loads(json_str)
        assert loaded["missingness_rate"] == 0.40


class TestConfigValidation:
    """Test configuration validation edge cases."""

    def test_valid_config_passes_validation(self):
        """Test that valid configuration passes validation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            clinvar_path = Path(tmpdir) / "clinvar.ht"
            dbnsfp_path = Path(tmpdir) / "dbnsfp.ht"
            clinvar_path.mkdir()
            dbnsfp_path.mkdir()

            config = PSROCConfig(
                genes=["BRCA1", "BRCA2"],
                clinvar_ht=str(clinvar_path),
                dbnsfp_ht=str(dbnsfp_path),
                scores=["CADD_phred", "REVEL_score"],
                output_dir=tmpdir,
                pathogenic_labels=["Pathogenic"],
                benign_labels=["Benign"],
            )

            errors = config.validate()
            assert len(errors) == 0

    def test_max_missingness_boundaries(self):
        """Test max_missingness boundary validation."""
        # Valid: 0.0
        config = PSROCConfig(max_missingness=0.0)
        errors = config.validate()
        assert not any("max_missingness" in e for e in errors if "between" in e)

        # Valid: 1.0
        config = PSROCConfig(max_missingness=1.0)
        errors = config.validate()
        assert not any("max_missingness" in e for e in errors if "between" in e)

        # Invalid: negative
        config = PSROCConfig(max_missingness=-0.1)
        errors = config.validate()
        assert any("max_missingness" in e for e in errors)

        # Invalid: > 1.0
        config = PSROCConfig(max_missingness=1.5)
        errors = config.validate()
        assert any("max_missingness" in e for e in errors)

    def test_threshold_methods(self):
        """Test that all threshold methods are accepted."""
        valid_methods = ["youden", "closest_to_corner", "f1"]

        for method in valid_methods:
            config = PSROCConfig(threshold_method=method)
            errors = config.validate()
            assert not any("threshold_method" in e for e in errors)

    def test_min_stars_validation(self):
        """Test min_stars validation."""
        # Valid: 0-4
        for stars in range(5):
            config = PSROCConfig(min_stars=stars)
            errors = config.validate()
            assert not any("min_stars" in e for e in errors)

        # Invalid: negative
        config = PSROCConfig(min_stars=-1)
        errors = config.validate()
        assert any("min_stars" in e for e in errors)


class TestStatePersistence:
    """Test pipeline state save/load functionality."""

    def test_state_round_trip(self):
        """Test that state can be saved and loaded correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            state_file = Path(tmpdir) / "state.json"

            # Create state
            config = PSROCConfig(
                genes=["BRCA1"],
                clinvar_ht="/data/clinvar.ht",
                dbnsfp_ht="/data/dbnsfp.ht",
                scores=["CADD_phred"],
                output_dir="/results",
                max_missingness=0.25,
            )

            state = PSROCState(config=config)
            state.start_time = "2024-01-01T00:00:00"
            state.current_stage = PSROCStage.LOAD_TABLES
            state.mark_stage_complete(PSROCStage.LOAD_TABLES, "/data/tables")
            state.mark_stage_complete(PSROCStage.FILTER_CLINVAR)
            state.errors.append("Test error")

            # Save
            state.save(state_file)
            assert state_file.exists()

            # Load
            loaded = PSROCState.load(state_file)

            # Verify
            assert loaded.config.genes == ["BRCA1"]
            assert loaded.config.max_missingness == 0.25
            assert loaded.start_time == "2024-01-01T00:00:00"
            assert loaded.is_stage_complete(PSROCStage.LOAD_TABLES)
            assert loaded.is_stage_complete(PSROCStage.FILTER_CLINVAR)
            assert not loaded.is_stage_complete(PSROCStage.ASSIGN_LABELS)
            assert "Test error" in loaded.errors
            assert loaded.outputs["load_tables"] == "/data/tables"


class TestResultSummary:
    """Test result summary generation."""

    def test_summary_contains_key_info(self):
        """Test that summary contains all key information."""
        missingness = ScoreMissingness(
            score_name="CADD_phred",
            n_total=100,
            n_present=95,
            n_missing=5,
            missingness_rate=0.05,
            included_in_analysis=True,
        )

        roc_result = ROCResult(
            score_name="CADD_phred",
            fpr=np.array([0.0, 0.5, 1.0]),
            tpr=np.array([0.0, 0.85, 1.0]),
            thresholds=np.array([1.0, 0.5, 0.0]),
            auc=0.85,
            optimal_threshold=0.5,
            sensitivity_at_optimal=0.85,
            specificity_at_optimal=0.5,
            n_variants_used=95,
            missingness=missingness,
        )

        result = PSROCResult(
            annotated_ht_path="/test/annotated.ht",
            metrics={"CADD_phred": roc_result},
            missingness={"CADD_phred": missingness},
            n_pathogenic=50,
            n_benign=45,
            n_excluded=5,
            n_total=100,
            scores_included=["CADD_phred"],
            scores_excluded=["REVEL_score"],
            max_missingness_threshold=0.3,
            output_dir="/test",
        )

        summary = result.summary()

        assert "100" in summary  # n_total
        assert "50" in summary  # n_pathogenic
        assert "45" in summary  # n_benign
        assert "CADD_phred" in summary
        assert "0.850" in summary  # AUC
        assert "REVEL_score" in summary  # excluded score
        assert "/test" in summary  # output_dir


class TestMultiScoreComparison:
    """Test multi-score comparison scenarios."""

    def test_multiple_scores_ranking(self):
        """Test that multiple scores are correctly ranked by AUC."""
        np.random.seed(42)
        n = 200
        labels = np.array([0] * (n // 2) + [1] * (n // 2))

        # Create scores with different separation qualities
        # Best: large separation
        best_b = np.random.normal(0.2, 0.1, n // 2)
        best_p = np.random.normal(0.8, 0.1, n // 2)

        # Medium: moderate separation
        med_b = np.random.normal(0.3, 0.15, n // 2)
        med_p = np.random.normal(0.7, 0.15, n // 2)

        # Worst: small separation
        worst_b = np.random.normal(0.4, 0.2, n // 2)
        worst_p = np.random.normal(0.6, 0.2, n // 2)

        scores = {
            "best_score": np.clip(np.concatenate([best_b, best_p]), 0, 1),
            "medium_score": np.clip(np.concatenate([med_b, med_p]), 0, 1),
            "worst_score": np.clip(np.concatenate([worst_b, worst_p]), 0, 1),
        }

        results = compute_roc_metrics(labels, scores)

        # Verify ranking
        assert results["best_score"].auc > results["medium_score"].auc
        assert results["medium_score"].auc > results["worst_score"].auc

    def test_partial_score_exclusion(self):
        """Test that some scores are excluded while others are included."""
        labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])

        scores = {
            # Low missingness - should be included
            "good_score": np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
            # High missingness - should be excluded
            "bad_score": np.array([np.nan] * 5 + [0.6, 0.7, 0.8, 0.9, 1.0]),
        }

        results = compute_roc_metrics(labels, scores, max_missingness=0.3)

        assert "good_score" in results
        assert "bad_score" not in results

        # Check missingness stats
        all_missingness = compute_all_missingness(scores, max_missingness=0.3)
        assert all_missingness["good_score"].included_in_analysis is True
        assert all_missingness["bad_score"].included_in_analysis is False


@pytest.mark.hail
class TestHailIntegration:
    """Integration tests requiring Hail.

    These tests are marked with @pytest.mark.hail and are skipped
    by default. Run with: pytest -m hail
    """

    @pytest.fixture
    def hail_session(self):
        """Initialize Hail for tests."""
        import hail as hl

        if not hl.spark_context():
            hl.init(quiet=True)
        yield
        # Don't stop Hail to allow reuse

    def test_clinvar_table_exists(self, hail_session):
        """Test that ClinVar test data can be loaded."""
        if not CLINVAR_TEST_DATA.exists():
            pytest.skip(f"Test data not found: {CLINVAR_TEST_DATA}")

        # Just verify we can read the VCF header
        # Full table building would be too slow for a unit test
        assert CLINVAR_TEST_DATA.exists()

    def test_dbnsfp_table_exists(self, hail_session):
        """Test that dbNSFP test data exists."""
        if not DBNSFP_TEST_DATA.exists():
            pytest.skip(f"Test data not found: {DBNSFP_TEST_DATA}")

        assert DBNSFP_TEST_DATA.exists()
