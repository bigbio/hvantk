"""Tests for ancestry HTML report generation.

Tests the report generation functions in hvantk.ancestry.report without
requiring Hail or actual genetic data. Uses mock objects to test report logic.
"""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple
from unittest.mock import MagicMock, patch

from hvantk.ancestry.constants import (
    PREDICTED_ANCESTRY_COL,
    ANCESTRY_PROB_COL,
    SOURCE_COL,
    KNOWN_ANCESTRY_COL,
)


@pytest.fixture
def sample_predictions_df():
    """Create sample predictions DataFrame for testing reports."""
    np.random.seed(42)
    data = []

    # Reference samples (30 each from EUR, AFR, EAS)
    for pop in ["EUR", "AFR", "EAS"]:
        for i in range(30):
            data.append(
                {
                    "s": f"ref_{pop}_{i}",
                    "PC1": np.random.normal(0, 1),
                    "PC2": np.random.normal(0, 1),
                    "PC3": np.random.normal(0, 1),
                    PREDICTED_ANCESTRY_COL: pop,
                    ANCESTRY_PROB_COL: np.random.uniform(0.85, 0.99),
                    KNOWN_ANCESTRY_COL: pop,
                    SOURCE_COL: "reference",
                }
            )

    # Query samples
    for i in range(50):
        assigned_pop = np.random.choice(["EUR", "AFR", "EAS"])
        prob = np.random.uniform(0.5, 0.95)
        predicted = assigned_pop if prob >= 0.75 else "unassigned"

        data.append(
            {
                "s": f"query_{i}",
                "PC1": np.random.normal(0, 2),
                "PC2": np.random.normal(0, 2),
                "PC3": np.random.normal(0, 1),
                PREDICTED_ANCESTRY_COL: predicted,
                ANCESTRY_PROB_COL: prob,
                KNOWN_ANCESTRY_COL: None,
                SOURCE_COL: "query",
            }
        )

    return pd.DataFrame(data)


@pytest.fixture
def mock_config():
    """Create mock pipeline configuration."""
    config = MagicMock()
    config.to_dict.return_value = {
        "min_af": 0.01,
        "max_af": 0.99,
        "min_call_rate": 0.98,
        "n_pcs": 20,
        "n_pcs_classify": 10,
        "min_prob": 0.75,
    }
    return config


@pytest.fixture
def mock_classification_result():
    """Create mock classification result."""
    np.random.seed(42)
    classes = ["EUR", "AFR", "EAS"]
    n_samples = 90

    y_true = np.array(classes * 30)
    y_pred = y_true.copy()
    # Simulate some mistakes
    mistakes = np.random.choice(n_samples, size=9, replace=False)
    for idx in mistakes:
        wrong_classes = [c for c in classes if c != y_true[idx]]
        y_pred[idx] = np.random.choice(wrong_classes)

    result = MagicMock()
    result.confusion_matrix = np.array(
        [
            [27, 2, 1],
            [1, 28, 1],
            [0, 2, 28],
        ]
    )
    result.confusion_matrix_labels = (y_true, y_pred)
    result.classes = classes
    result.validation_metrics = {"accuracy": 0.9}
    return result


@pytest.fixture
def mock_ancestry_result(
    sample_predictions_df, mock_config, mock_classification_result
):
    """Create mock AncestryInferenceResult."""
    result = MagicMock()
    result.get_predictions_df.return_value = sample_predictions_df
    result.config = mock_config
    result.classification_result = mock_classification_result
    result.eigenvalues = [10 * (0.7**i) for i in range(20)]
    result.pipeline_stats = {
        "n_query_samples": 50,
        "n_reference_samples": 90,
        "n_populations": 3,
        "n_shared_variants": 1000,
    }

    def get_accuracy():
        if mock_classification_result.validation_metrics:
            return mock_classification_result.validation_metrics.get("accuracy")
        return None

    result.get_accuracy = get_accuracy
    return result


class TestCreateSummaryCard:
    """Tests for _create_summary_card function."""

    def test_create_summary_card_basic(self):
        """Test creating a basic summary card."""
        from hvantk.ancestry.report import _create_summary_card

        html = _create_summary_card(100, "Total Samples")
        assert "100" in html
        assert "Total Samples" in html
        assert 'class="summary-card"' in html

    def test_create_summary_card_formatted_number(self):
        """Test summary card with formatted number."""
        from hvantk.ancestry.report import _create_summary_card

        html = _create_summary_card("1,234", "Large Count")
        assert "1,234" in html

    def test_create_summary_card_percentage(self):
        """Test summary card with percentage."""
        from hvantk.ancestry.report import _create_summary_card

        html = _create_summary_card("95.5%", "Accuracy")
        assert "95.5%" in html


class TestCreateAncestryTable:
    """Tests for _create_ancestry_table function."""

    def test_create_ancestry_table(self, sample_predictions_df):
        """Test creating ancestry distribution table."""
        from hvantk.ancestry.report import _create_ancestry_table

        html = _create_ancestry_table(sample_predictions_df)

        assert "<table>" in html
        assert "</table>" in html
        assert "Population" in html or "Code" in html

    def test_ancestry_table_has_populations(self, sample_predictions_df):
        """Test that table includes population names."""
        from hvantk.ancestry.report import _create_ancestry_table

        html = _create_ancestry_table(sample_predictions_df)

        # Should contain at least one population
        assert any(pop in html for pop in ["EUR", "AFR", "EAS", "unassigned"])


class TestCreatePredictionsTable:
    """Tests for _create_predictions_table function."""

    def test_create_predictions_table(self, sample_predictions_df):
        """Test creating sample predictions table."""
        from hvantk.ancestry.report import _create_predictions_table

        html = _create_predictions_table(sample_predictions_df)

        assert "<table>" in html
        assert "</table>" in html

    def test_predictions_table_max_rows(self, sample_predictions_df):
        """Test predictions table respects max_rows."""
        from hvantk.ancestry.report import _create_predictions_table

        html = _create_predictions_table(sample_predictions_df, max_rows=10)

        # Should mention limited rows
        assert "10" in html or "Showing" in html.lower()

    def test_predictions_table_includes_sample_ids(self, sample_predictions_df):
        """Test that table includes sample IDs."""
        from hvantk.ancestry.report import _create_predictions_table

        html = _create_predictions_table(sample_predictions_df, max_rows=5)

        # Should have sample ID in header or rows
        assert "s" in html or "sample" in html.lower()


class TestCreateConfigTable:
    """Tests for _create_config_table function."""

    def test_create_config_table(self):
        """Test creating configuration table."""
        from hvantk.ancestry.report import _create_config_table

        config_dict = {
            "min_af": 0.01,
            "max_af": 0.99,
            "n_pcs": 20,
        }
        html = _create_config_table(config_dict)

        assert "<table" in html
        assert "min_af" in html
        assert "0.01" in html

    def test_config_table_all_params(self):
        """Test config table includes all parameters."""
        from hvantk.ancestry.report import _create_config_table

        config_dict = {
            "param1": "value1",
            "param2": 123,
            "param3": True,
        }
        html = _create_config_table(config_dict)

        for key in config_dict:
            assert key in html


class TestCreateValidationSection:
    """Tests for _create_validation_section function."""

    def test_validation_section_with_metrics(self, mock_ancestry_result):
        """Test validation section when metrics are present."""
        from hvantk.ancestry.report import _create_validation_section

        html = _create_validation_section(mock_ancestry_result)

        # Should show accuracy
        assert "90" in html or "0.9" in html or "Accuracy" in html.lower()

    def test_validation_section_skipped(self):
        """Test validation section when validation was skipped."""
        from hvantk.ancestry.report import _create_validation_section

        result = MagicMock()
        result.get_accuracy.return_value = None
        result.classification_result.confusion_matrix = None

        html = _create_validation_section(result)

        # Should show message about skipped validation
        assert "skip" in html.lower() or "enable" in html.lower()


class TestGenerateAncestryReport:
    """Tests for generate_ancestry_report function."""

    def test_generate_report_creates_file(self, mock_ancestry_result, tmp_path):
        """Test that report generation creates HTML file."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "report.html"
        result_path = generate_ancestry_report(mock_ancestry_result, output_path)

        assert result_path == output_path
        assert output_path.exists()

    def test_report_contains_html_structure(self, mock_ancestry_result, tmp_path):
        """Test that generated report has proper HTML structure."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "report.html"
        generate_ancestry_report(mock_ancestry_result, output_path)

        content = output_path.read_text()
        assert "<!DOCTYPE html>" in content
        assert "<html" in content
        assert "</html>" in content
        assert "<head>" in content
        assert "<body>" in content

    def test_report_contains_title(self, mock_ancestry_result, tmp_path):
        """Test that report contains custom title."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "report.html"
        generate_ancestry_report(
            mock_ancestry_result,
            output_path,
            title="My Custom Report",
        )

        content = output_path.read_text()
        assert "My Custom Report" in content

    def test_report_contains_summary(self, mock_ancestry_result, tmp_path):
        """Test that report contains summary section."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "report.html"
        generate_ancestry_report(mock_ancestry_result, output_path)

        content = output_path.read_text()
        assert "Summary" in content
        assert "50" in content  # n_query_samples

    def test_report_contains_ancestry_section(self, mock_ancestry_result, tmp_path):
        """Test that report contains ancestry distribution section."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "report.html"
        generate_ancestry_report(mock_ancestry_result, output_path)

        content = output_path.read_text()
        assert "Ancestry" in content

    def test_report_contains_pca_section(self, mock_ancestry_result, tmp_path):
        """Test that report contains PCA visualization section."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "report.html"
        generate_ancestry_report(mock_ancestry_result, output_path)

        content = output_path.read_text()
        assert "PCA" in content

    def test_report_contains_embedded_images(self, mock_ancestry_result, tmp_path):
        """Test that report contains base64 encoded images."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "report.html"
        generate_ancestry_report(mock_ancestry_result, output_path)

        content = output_path.read_text()
        # Should have embedded images
        assert "data:image/png;base64" in content

    def test_report_contains_config_section(self, mock_ancestry_result, tmp_path):
        """Test that report contains configuration section."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "report.html"
        generate_ancestry_report(mock_ancestry_result, output_path)

        content = output_path.read_text()
        assert "Configuration" in content or "Config" in content
        assert "min_af" in content

    def test_report_creates_parent_directory(self, mock_ancestry_result, tmp_path):
        """Test that report creates parent directories if needed."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "subdir" / "nested" / "report.html"
        generate_ancestry_report(mock_ancestry_result, output_path)

        assert output_path.exists()

    def test_report_returns_path(self, mock_ancestry_result, tmp_path):
        """Test that function returns Path object."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "report.html"
        result = generate_ancestry_report(mock_ancestry_result, output_path)

        assert isinstance(result, Path)

    def test_report_with_string_path(self, mock_ancestry_result, tmp_path):
        """Test that function accepts string path."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = str(tmp_path / "report.html")
        generate_ancestry_report(mock_ancestry_result, output_path)

        assert Path(output_path).exists()


class TestReportStyling:
    """Tests for report CSS and styling."""

    def test_report_has_css(self, mock_ancestry_result, tmp_path):
        """Test that report includes CSS styles."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "report.html"
        generate_ancestry_report(mock_ancestry_result, output_path)

        content = output_path.read_text()
        assert "<style>" in content
        assert "</style>" in content

    def test_report_has_ancestry_badges(self, mock_ancestry_result, tmp_path):
        """Test that report has ancestry badge styling."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "report.html"
        generate_ancestry_report(mock_ancestry_result, output_path)

        content = output_path.read_text()
        assert "ancestry-badge" in content or "badge" in content.lower()

    def test_report_has_responsive_design(self, mock_ancestry_result, tmp_path):
        """Test that report has responsive design elements."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "report.html"
        generate_ancestry_report(mock_ancestry_result, output_path)

        content = output_path.read_text()
        assert "viewport" in content


class TestReportFooter:
    """Tests for report footer."""

    def test_report_has_footer(self, mock_ancestry_result, tmp_path):
        """Test that report has footer section."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "report.html"
        generate_ancestry_report(mock_ancestry_result, output_path)

        content = output_path.read_text()
        assert "footer" in content.lower() or "hvantk" in content

    def test_report_has_timestamp(self, mock_ancestry_result, tmp_path):
        """Test that report includes generation timestamp."""
        from hvantk.ancestry.report import generate_ancestry_report

        output_path = tmp_path / "report.html"
        generate_ancestry_report(mock_ancestry_result, output_path)

        content = output_path.read_text()
        # Should have date-like pattern
        assert "Generated" in content or "20" in content  # Year prefix
