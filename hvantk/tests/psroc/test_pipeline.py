"""
Tests for PSROC pipeline module.

This module tests the pipeline orchestration including configuration validation,
label assignment, and end-to-end pipeline execution.
"""

import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pytest

from hvantk.algorithms.psroc.pipeline import (
    PSROCConfig,
    PSROCState,
    PSROCResult,
    PSROCPipeline,
    PSROCStage,
    PATHOGENIC_LABELS,
    BENIGN_LABELS,
    CLNREVSTAT_STAR_MAP,
    SCORE_DIRECTIONALITY,
    parse_variant_list,
)
from hvantk.algorithms.psroc.roc import ROCResult, ScoreMissingness


class TestPSROCConfig:
    """Test PSROCConfig dataclass and validation."""

    def test_config_creation_with_genes(self):
        """Test creating a config with gene list."""
        config = PSROCConfig(
            genes=["BRCA1", "BRCA2"],
            clinvar_ht="/data/clinvar.ht",
            dbnsfp_ht="/data/dbnsfp.ht",
            scores=["CADD_phred", "REVEL_score"],
            output_dir="/results/psroc",
        )

        assert config.genes == ["BRCA1", "BRCA2"]
        assert config.clinvar_ht == "/data/clinvar.ht"
        assert config.scores == ["CADD_phred", "REVEL_score"]
        assert config.max_missingness == 0.3  # default
        assert config.min_stars == 1  # default

    def test_config_default_values(self):
        """Test default configuration values."""
        config = PSROCConfig()

        assert config.reference_genome == "GRCh38"
        assert config.min_stars == 1
        assert config.max_missingness == 0.3
        assert config.threshold_method == "youden"
        assert config.export_tsv is False
        assert config.overwrite is False
        assert config.generate_plots is True

    def test_validation_missing_variant_source(self):
        """Test validation fails when no variant source provided."""
        config = PSROCConfig(
            clinvar_ht="/data/clinvar.ht",
            dbnsfp_ht="/data/dbnsfp.ht",
            scores=["CADD_phred"],
            output_dir="/results",
        )

        errors = config.validate()
        assert any("--genes" in e for e in errors)

    def test_validation_multiple_variant_sources(self):
        """Test validation fails with multiple variant sources."""
        config = PSROCConfig(
            genes=["BRCA1"],
            variants_path="/data/variants.txt",
            clinvar_ht="/data/clinvar.ht",
            dbnsfp_ht="/data/dbnsfp.ht",
            scores=["CADD_phred"],
            output_dir="/results",
        )

        errors = config.validate()
        assert any("Cannot provide multiple variant sources" in e for e in errors)

    def test_config_with_gene_set_collection(self):
        """Test config accepts gene_set_collection as valid source."""
        collection = {
            "cardiac": {"BRCA1", "TP53"},
            "neuro": {"SCN1A", "SCN2A"},
        }
        config = PSROCConfig(
            gene_set_collection=collection,
            clinvar_ht="/data/clinvar.ht",
            dbnsfp_ht="/data/dbnsfp.ht",
            scores=["CADD_phred"],
            output_dir="/results",
        )

        assert config.gene_set_collection == collection
        # Should not produce source-related errors (table paths won't exist)
        errors = config.validate()
        assert not any("--genes" in e and "Must provide" in e for e in errors)

    def test_validation_gene_set_collection_with_genes_is_error(self):
        """Test validation rejects gene_set_collection + genes together."""
        config = PSROCConfig(
            genes=["BRCA1"],
            gene_set_collection={"cardiac": {"TP53"}},
            clinvar_ht="/data/clinvar.ht",
            dbnsfp_ht="/data/dbnsfp.ht",
            scores=["CADD_phred"],
            output_dir="/results",
        )

        errors = config.validate()
        assert any("Cannot provide multiple variant sources" in e for e in errors)

    def test_validation_gene_set_collection_empty_groups(self):
        """Test validation catches empty groups in gene_set_collection."""
        config = PSROCConfig(
            gene_set_collection={"cardiac": {"BRCA1"}, "empty": set()},
            clinvar_ht="/data/clinvar.ht",
            dbnsfp_ht="/data/dbnsfp.ht",
            scores=["CADD_phred"],
            output_dir="/results",
        )

        errors = config.validate()
        assert any("empty groups" in e for e in errors)

    def test_validation_missing_clinvar_path(self):
        """Test validation fails when clinvar_ht not provided."""
        config = PSROCConfig(
            genes=["BRCA1"],
            dbnsfp_ht="/data/dbnsfp.ht",
            scores=["CADD_phred"],
            output_dir="/results",
        )

        errors = config.validate()
        assert any("--clinvar-ht is required" in e for e in errors)

    def test_validation_missing_dbnsfp_path(self):
        """Test validation fails when dbnsfp_ht not provided."""
        config = PSROCConfig(
            genes=["BRCA1"],
            clinvar_ht="/data/clinvar.ht",
            scores=["CADD_phred"],
            output_dir="/results",
        )

        errors = config.validate()
        assert any("--dbnsfp-ht is required" in e for e in errors)

    def test_validation_missing_scores(self):
        """Test validation fails when no scores provided."""
        config = PSROCConfig(
            genes=["BRCA1"],
            clinvar_ht="/data/clinvar.ht",
            dbnsfp_ht="/data/dbnsfp.ht",
            output_dir="/results",
        )

        errors = config.validate()
        assert any("--scores" in e for e in errors)

    def test_validation_missing_output_dir(self):
        """Test validation fails when output_dir not provided."""
        config = PSROCConfig(
            genes=["BRCA1"],
            clinvar_ht="/data/clinvar.ht",
            dbnsfp_ht="/data/dbnsfp.ht",
            scores=["CADD_phred"],
        )

        errors = config.validate()
        assert any("--output-dir is required" in e for e in errors)

    def test_validation_invalid_max_missingness(self):
        """Test validation fails with invalid max_missingness."""
        config = PSROCConfig(
            genes=["BRCA1"],
            clinvar_ht="/data/clinvar.ht",
            dbnsfp_ht="/data/dbnsfp.ht",
            scores=["CADD_phred"],
            output_dir="/results",
            max_missingness=1.5,  # Invalid: > 1.0
        )

        errors = config.validate()
        assert any("max_missingness" in e for e in errors)

    def test_validation_invalid_threshold_method(self):
        """Test validation fails with invalid threshold method."""
        config = PSROCConfig(
            genes=["BRCA1"],
            clinvar_ht="/data/clinvar.ht",
            dbnsfp_ht="/data/dbnsfp.ht",
            scores=["CADD_phred"],
            output_dir="/results",
            threshold_method="invalid_method",
        )

        errors = config.validate()
        assert any("Invalid threshold_method" in e for e in errors)

    def test_get_gene_set_from_genes(self):
        """Test get_gene_set returns set from genes list."""
        config = PSROCConfig(genes=["BRCA1", "BRCA2", "TP53"])

        gene_set = config.get_gene_set()

        assert gene_set == {"BRCA1", "BRCA2", "TP53"}

    def test_get_gene_set_with_variants_path(self):
        """Test get_gene_set returns None when using variants_path."""
        config = PSROCConfig(variants_path="/data/variants.txt")

        gene_set = config.get_gene_set()

        assert gene_set is None

    def test_get_gene_set_from_file(self):
        """Test get_gene_set loads genes from file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("BRCA1\n")
            f.write("BRCA2\n")
            f.write("# comment\n")
            f.write("TP53\n")
            genes_file = f.name

        try:
            config = PSROCConfig(genes_file=genes_file)
            gene_set = config.get_gene_set()

            assert "BRCA1" in gene_set
            assert "BRCA2" in gene_set
            assert "TP53" in gene_set
            assert len(gene_set) == 3
        finally:
            Path(genes_file).unlink()


class TestPSROCState:
    """Test PSROCState dataclass."""

    def test_state_creation(self):
        """Test creating a state instance."""
        config = PSROCConfig(
            genes=["BRCA1"],
            clinvar_ht="/data/clinvar.ht",
            dbnsfp_ht="/data/dbnsfp.ht",
            scores=["CADD_phred"],
            output_dir="/results",
        )

        state = PSROCState(config=config)

        assert state.config == config
        assert state.current_stage is None
        assert state.completed_stages == []
        assert state.outputs == {}
        assert state.errors == []

    def test_mark_stage_complete(self):
        """Test marking a stage as complete."""
        config = PSROCConfig()
        state = PSROCState(config=config)

        state.mark_stage_complete(PSROCStage.LOAD_TABLES, "/data/tables")

        assert "load_tables" in state.completed_stages
        assert state.outputs["load_tables"] == "/data/tables"

    def test_is_stage_complete(self):
        """Test checking stage completion status."""
        config = PSROCConfig()
        state = PSROCState(config=config)

        assert not state.is_stage_complete(PSROCStage.LOAD_TABLES)

        state.mark_stage_complete(PSROCStage.LOAD_TABLES)

        assert state.is_stage_complete(PSROCStage.LOAD_TABLES)

    def test_save_and_load(self):
        """Test state serialization and deserialization."""
        config = PSROCConfig(
            genes=["BRCA1"],
            clinvar_ht="/data/clinvar.ht",
            dbnsfp_ht="/data/dbnsfp.ht",
            scores=["CADD_phred"],
            output_dir="/results",
        )
        state = PSROCState(config=config)
        state.mark_stage_complete(PSROCStage.LOAD_TABLES, "/data/tables")
        state.start_time = "2024-01-01T00:00:00"

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            state_file = Path(f.name)

        try:
            state.save(state_file)
            loaded_state = PSROCState.load(state_file)

            assert loaded_state.config.genes == ["BRCA1"]
            assert loaded_state.completed_stages == ["load_tables"]
            assert loaded_state.outputs["load_tables"] == "/data/tables"
            assert loaded_state.start_time == "2024-01-01T00:00:00"
        finally:
            state_file.unlink()


class TestPSROCResult:
    """Test PSROCResult dataclass."""

    def _create_mock_result(self) -> PSROCResult:
        """Create a mock result for testing."""
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

        return PSROCResult(
            annotated_ht_path="/results/annotated.ht",
            metrics={"CADD_phred": roc_result},
            missingness={"CADD_phred": missingness},
            n_pathogenic=50,
            n_benign=45,
            n_excluded=5,
            n_total=100,
            scores_included=["CADD_phred"],
            scores_excluded=[],
            max_missingness_threshold=0.3,
            output_dir="/results",
        )

    def test_result_creation(self):
        """Test creating a result instance."""
        result = self._create_mock_result()

        assert result.n_pathogenic == 50
        assert result.n_benign == 45
        assert result.n_excluded == 5
        assert result.n_total == 100
        assert "CADD_phred" in result.metrics
        assert result.metrics["CADD_phred"].auc == 0.85

    def test_result_to_dict(self):
        """Test result serialization to dict."""
        result = self._create_mock_result()

        d = result.to_dict()

        assert d["n_pathogenic"] == 50
        assert d["n_benign"] == 45
        assert d["n_total"] == 100
        assert "CADD_phred" in d["metrics"]
        assert d["metrics"]["CADD_phred"]["auc"] == 0.85
        assert d["scores_included"] == ["CADD_phred"]

    def test_result_summary(self):
        """Test result summary generation."""
        result = self._create_mock_result()

        summary = result.summary()

        assert "Variants analyzed: 100" in summary
        assert "Pathogenic (P): 50" in summary
        assert "Benign (B): 45" in summary
        assert "CADD_phred" in summary
        assert "AUC=0.850" in summary


class TestPSROCStage:
    """Test PSROCStage enum."""

    def test_stage_values(self):
        """Test stage enum values."""
        assert PSROCStage.LOAD_TABLES.value == "load_tables"
        assert PSROCStage.FILTER_CLINVAR.value == "filter_clinvar"
        assert PSROCStage.ASSIGN_LABELS.value == "assign_labels"
        assert PSROCStage.ANNOTATE_SCORES.value == "annotate_scores"
        assert PSROCStage.COMPUTE_MISSINGNESS.value == "compute_missingness"
        assert PSROCStage.COMPUTE_ROC.value == "compute_roc"
        assert PSROCStage.GENERATE_OUTPUTS.value == "generate_outputs"


class TestConstants:
    """Test pipeline constants (labels, star map)."""

    def test_pathogenic_labels(self):
        """Test pathogenic label list."""
        assert "Pathogenic" in PATHOGENIC_LABELS
        assert "Likely_pathogenic" in PATHOGENIC_LABELS
        assert "Pathogenic/Likely_pathogenic" in PATHOGENIC_LABELS
        assert len(PATHOGENIC_LABELS) == 3

    def test_benign_labels(self):
        """Test benign label list."""
        assert "Benign" in BENIGN_LABELS
        assert "Likely_benign" in BENIGN_LABELS
        assert "Benign/Likely_benign" in BENIGN_LABELS
        assert len(BENIGN_LABELS) == 3

    def test_clnrevstat_star_map_known_values(self):
        """Test CLNREVSTAT star map assigns correct star counts."""
        assert CLNREVSTAT_STAR_MAP["practice_guideline"] == 4
        assert CLNREVSTAT_STAR_MAP["reviewed_by_expert_panel"] == 3
        assert (
            CLNREVSTAT_STAR_MAP["criteria_provided,_multiple_submitters,_no_conflicts"]
            == 2
        )
        assert CLNREVSTAT_STAR_MAP["criteria_provided,_single_submitter"] == 1
        assert CLNREVSTAT_STAR_MAP["no_assertion_criteria_provided"] == 0

    def test_clnrevstat_star_map_values_in_range(self):
        """Test all star values are valid (0-4)."""
        for status, stars in CLNREVSTAT_STAR_MAP.items():
            assert 0 <= stars <= 4, f"{status} has invalid star count {stars}"

    def test_clnrevstat_star_map_conflicting_classifications(self):
        """Test both ClinVar conflicting status variants are mapped."""
        # ClinVar has used both naming conventions across versions
        assert "criteria_provided,_conflicting_classifications" in CLNREVSTAT_STAR_MAP
        assert "criteria_provided,_conflicting_interpretations" in CLNREVSTAT_STAR_MAP
        assert (
            CLNREVSTAT_STAR_MAP["criteria_provided,_conflicting_classifications"]
            == CLNREVSTAT_STAR_MAP["criteria_provided,_conflicting_interpretations"]
        )

    def test_score_directionality_higher_is_pathogenic(self):
        """Test higher-is-pathogenic scores are marked True."""
        assert SCORE_DIRECTIONALITY["CADD_phred"] is True
        assert SCORE_DIRECTIONALITY["REVEL_score"] is True
        assert SCORE_DIRECTIONALITY["MetaLR_score"] is True

    def test_score_directionality_lower_is_pathogenic(self):
        """Test lower-is-pathogenic scores are marked False."""
        assert SCORE_DIRECTIONALITY["SIFT_score"] is False
        assert SCORE_DIRECTIONALITY["SIFT4G_score"] is False
        assert SCORE_DIRECTIONALITY["PROVEAN_score"] is False
        assert SCORE_DIRECTIONALITY["FATHMM_score"] is False
        assert SCORE_DIRECTIONALITY["LRT_score"] is False

    def test_score_directionality_values_are_bool(self):
        """Test all direction values are bool."""
        for score, direction in SCORE_DIRECTIONALITY.items():
            assert isinstance(direction, bool), f"{score} has non-bool direction"

    def test_config_score_directions_override(self):
        """Test user-provided score_directions override defaults."""
        config = PSROCConfig(
            genes=["BRCA1"],
            clinvar_ht="/data/clinvar.ht",
            dbnsfp_ht="/data/dbnsfp.ht",
            scores=["CADD_phred", "SIFT_score"],
            output_dir="/results",
            score_directions={"SIFT_score": True},  # override
        )
        assert config.score_directions == {"SIFT_score": True}


class TestPSROCPipelineValidation:
    """Test PSROCPipeline initialization and validation."""

    def test_pipeline_invalid_config_raises_error(self):
        """Test that invalid config raises ValueError on init."""
        config = PSROCConfig(
            # Missing required fields
            genes=["BRCA1"],
            output_dir="/results",
        )

        with pytest.raises(ValueError, match="Configuration validation failed"):
            PSROCPipeline(config)

    @patch("hvantk.algorithms.psroc.pipeline.hl")
    @patch("hvantk.algorithms.psroc.pipeline.Path.exists")
    def test_pipeline_setup_output_paths(self, mock_exists, mock_hl):
        """Test that output paths are set up correctly."""
        mock_exists.return_value = True
        mock_hl.current_backend.return_value = Mock()

        with tempfile.TemporaryDirectory() as tmpdir:
            config = PSROCConfig(
                genes=["BRCA1"],
                clinvar_ht=f"{tmpdir}/clinvar.ht",
                dbnsfp_ht=f"{tmpdir}/dbnsfp.ht",
                scores=["CADD_phred"],
                output_dir=tmpdir,
                output_prefix="test",
            )

            # Create mock table directories
            Path(f"{tmpdir}/clinvar.ht").mkdir()
            Path(f"{tmpdir}/dbnsfp.ht").mkdir()

            pipeline = PSROCPipeline(config)

            assert "annotated_ht" in pipeline.paths
            assert "metrics_json" in pipeline.paths
            assert "roc_curves_png" in pipeline.paths
            assert pipeline.paths["annotated_ht"].endswith("test_annotated.ht")

    @patch("hvantk.algorithms.psroc.pipeline.hl")
    @patch("hvantk.algorithms.psroc.pipeline.Path.exists")
    def test_run_raises_when_gene_set_collection(self, mock_exists, mock_hl):
        """Test run() raises ValueError when gene_set_collection is active."""
        mock_exists.return_value = True
        mock_hl.current_backend.return_value = Mock()

        with tempfile.TemporaryDirectory() as tmpdir:
            config = PSROCConfig(
                gene_set_collection={"cardiac": {"BRCA1", "TP53"}},
                clinvar_ht=f"{tmpdir}/clinvar.ht",
                dbnsfp_ht=f"{tmpdir}/dbnsfp.ht",
                scores=["CADD_phred"],
                output_dir=tmpdir,
            )

            Path(f"{tmpdir}/clinvar.ht").mkdir()
            Path(f"{tmpdir}/dbnsfp.ht").mkdir()

            pipeline = PSROCPipeline(config)

            with pytest.raises(ValueError, match="run_collection"):
                pipeline.run()

    @patch("hvantk.algorithms.psroc.pipeline.hl")
    @patch("hvantk.algorithms.psroc.pipeline.Path.exists")
    def test_run_collection_raises_without_collection(self, mock_exists, mock_hl):
        """Test run_collection() raises when gene_set_collection not set."""
        mock_exists.return_value = True
        mock_hl.current_backend.return_value = Mock()

        with tempfile.TemporaryDirectory() as tmpdir:
            config = PSROCConfig(
                genes=["BRCA1"],
                clinvar_ht=f"{tmpdir}/clinvar.ht",
                dbnsfp_ht=f"{tmpdir}/dbnsfp.ht",
                scores=["CADD_phred"],
                output_dir=tmpdir,
            )

            Path(f"{tmpdir}/clinvar.ht").mkdir()
            Path(f"{tmpdir}/dbnsfp.ht").mkdir()

            pipeline = PSROCPipeline(config)

            with pytest.raises(
                ValueError, match="gene_set_collection is not configured"
            ):
                pipeline.run_collection()


class TestVariantFileParsing:
    """Test variant file parsing in the pipeline."""

    def test_parse_variant_file_format(self):
        """Test parsing variant file format chr:pos:ref:alt using the real parser."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("chr1:12345:A:T\n")
            f.write("chr2:67890:G:C\n")
            f.write("# comment line\n")
            f.write("chr3:11111:T:A\n")
            variants_file = f.name

        try:
            variants = parse_variant_list(variants_file)
            assert len(variants) == 3
            assert variants[0] == ("chr1", 12345, "A", "T")
            assert variants[1] == ("chr2", 67890, "G", "C")
            assert variants[2] == ("chr3", 11111, "T", "A")
        finally:
            Path(variants_file).unlink()

    def test_parse_variant_file_without_chr_prefix(self):
        """Test parsing variant file with variants without chr prefix using the real parser."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("1:12345:A:T\n")  # No chr prefix
            f.write("2:67890:G:C\n")
            variants_file = f.name

        try:
            variants = parse_variant_list(variants_file)
            assert len(variants) == 2
            assert variants[0][0] == "chr1"  # chr prefix added
            assert variants[1][0] == "chr2"
            assert variants[0][1] == 12345
            assert variants[1][1] == 67890
        finally:
            Path(variants_file).unlink()


class TestShowPlan:
    """Test show_plan method output."""

    @patch("hvantk.algorithms.psroc.pipeline.hl")
    @patch("hvantk.algorithms.psroc.pipeline.Path.exists")
    def test_show_plan_with_genes(self, mock_exists, mock_hl, capsys):
        """Test show_plan output with gene input."""
        mock_exists.return_value = True
        mock_hl.current_backend.return_value = Mock()

        with tempfile.TemporaryDirectory() as tmpdir:
            config = PSROCConfig(
                genes=["BRCA1", "BRCA2"],
                clinvar_ht=f"{tmpdir}/clinvar.ht",
                dbnsfp_ht=f"{tmpdir}/dbnsfp.ht",
                scores=["CADD_phred", "REVEL_score"],
                output_dir=tmpdir,
            )

            Path(f"{tmpdir}/clinvar.ht").mkdir()
            Path(f"{tmpdir}/dbnsfp.ht").mkdir()

            pipeline = PSROCPipeline(config)
            pipeline.show_plan()

            captured = capsys.readouterr()

            assert "PSROC PIPELINE EXECUTION PLAN" in captured.out
            assert "BRCA1" in captured.out
            assert "BRCA2" in captured.out
            assert "CADD_phred" in captured.out
            assert "REVEL_score" in captured.out
            assert "Load Hail Tables" in captured.out
            assert "Compute ROC" in captured.out
