import json
import os
from dataclasses import dataclass
from unittest.mock import MagicMock, patch

import pytest
import yaml

TESTDATA_DIR = os.path.join(os.path.dirname(__file__), "testdata")
CONFIG_PATH = os.path.join(TESTDATA_DIR, "alphagenome_config.yaml")


class TestLoadConfig:
    def test_load_valid_config(self):
        from hvantk.data.alphagenome_streamer import load_config

        config = load_config(CONFIG_PATH)
        assert config["api"]["key"] == "test-api-key-123"
        assert config["api"]["max_retries"] == 3
        assert config["ontology"]["terms"] == ["UBERON:0001157", "UBERON:0000955"]
        assert config["ontology"]["output_types"] == ["RNA_SEQ", "CHROMATIN"]
        assert config["intervals"]["adaptive"] is True

    def test_load_config_file_not_found(self):
        from hvantk.data.alphagenome_streamer import load_config

        with pytest.raises(FileNotFoundError):
            load_config("/nonexistent/path.yaml")

    def test_load_config_missing_api_section(self, tmp_path):
        from hvantk.data.alphagenome_streamer import load_config

        bad_config = tmp_path / "bad.yaml"
        bad_config.write_text(yaml.dump({"ontology": {"terms": []}}))
        with pytest.raises(ValueError, match="api"):
            load_config(str(bad_config))

    def test_load_config_missing_ontology_section(self, tmp_path):
        from hvantk.data.alphagenome_streamer import load_config

        bad_config = tmp_path / "bad.yaml"
        bad_config.write_text(yaml.dump({"api": {"key": "x"}, "intervals": {}}))
        with pytest.raises(ValueError, match="ontology"):
            load_config(str(bad_config))

    def test_resolve_api_key_from_env(self, tmp_path, monkeypatch):
        from hvantk.data.alphagenome_streamer import load_config

        cfg = {
            "api": {"key": None, "max_retries": 3, "retry_backoff": 2.0, "request_timeout": 120},
            "ontology": {"terms": ["UBERON:0001157"], "output_types": ["RNA_SEQ"]},
            "intervals": {"default_size": 1048576, "adaptive": True,
                          "adaptive_max_size": 1048576, "density_window": 50000},
        }
        cfg_path = tmp_path / "cfg.yaml"
        cfg_path.write_text(yaml.dump(cfg))
        monkeypatch.setenv("ALPHAGENOME_API_KEY", "env-key-456")
        config = load_config(str(cfg_path))
        assert config["api"]["key"] == "env-key-456"

    def test_resolve_api_key_missing_everywhere(self, tmp_path, monkeypatch):
        from hvantk.data.alphagenome_streamer import load_config

        cfg = {
            "api": {"key": None, "max_retries": 3, "retry_backoff": 2.0, "request_timeout": 120},
            "ontology": {"terms": ["UBERON:0001157"], "output_types": ["RNA_SEQ"]},
            "intervals": {"default_size": 1048576, "adaptive": True,
                          "adaptive_max_size": 1048576, "density_window": 50000},
        }
        cfg_path = tmp_path / "cfg.yaml"
        cfg_path.write_text(yaml.dump(cfg))
        monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
        with pytest.raises(ValueError, match="API key"):
            load_config(str(cfg_path))


@dataclass
class SimpleVariant:
    chrom: str
    pos: int
    ref: str
    alt: str


class TestComputeIntervals:
    """Test compute_intervals() as a pure function."""

    def _make_config(self, adaptive=True, default_size=1_048_576,
                     adaptive_max_size=1_048_576, density_window=50_000):
        return {
            "intervals": {
                "adaptive": adaptive,
                "default_size": default_size,
                "adaptive_max_size": adaptive_max_size,
                "density_window": density_window,
            }
        }

    def test_single_variant_fixed_mode(self):
        from hvantk.data.alphagenome_streamer import compute_intervals

        variants = [SimpleVariant("chr1", 500_000, "A", "T")]
        config = self._make_config(adaptive=False, default_size=100_000)
        result = compute_intervals(variants, config)
        assert len(result) == 1
        interval, grouped_variants = result[0]
        assert interval.chrom == "chr1"
        assert interval.start <= 500_000
        assert interval.end >= 500_000
        assert interval.end - interval.start == 100_000
        assert len(grouped_variants) == 1

    def test_two_close_variants_adaptive_grouped(self):
        from hvantk.data.alphagenome_streamer import compute_intervals

        variants = [
            SimpleVariant("chr1", 100_000, "A", "T"),
            SimpleVariant("chr1", 120_000, "G", "C"),
        ]
        config = self._make_config(adaptive=True, density_window=50_000)
        result = compute_intervals(variants, config)
        assert len(result) == 1
        assert len(result[0][1]) == 2

    def test_two_distant_variants_adaptive_separate(self):
        from hvantk.data.alphagenome_streamer import compute_intervals

        variants = [
            SimpleVariant("chr1", 100_000, "A", "T"),
            SimpleVariant("chr1", 2_000_000, "G", "C"),
        ]
        config = self._make_config(adaptive=True, density_window=50_000)
        result = compute_intervals(variants, config)
        assert len(result) == 2
        assert len(result[0][1]) == 1
        assert len(result[1][1]) == 1

    def test_different_chromosomes_always_separate(self):
        from hvantk.data.alphagenome_streamer import compute_intervals

        variants = [
            SimpleVariant("chr1", 100_000, "A", "T"),
            SimpleVariant("chr2", 100_010, "G", "C"),
        ]
        config = self._make_config(adaptive=True, density_window=50_000)
        result = compute_intervals(variants, config)
        assert len(result) == 2

    def test_large_cluster_split_at_max_size(self):
        from hvantk.data.alphagenome_streamer import compute_intervals

        variants = [SimpleVariant("chr1", i * 10_000, "A", "T") for i in range(200)]
        config = self._make_config(
            adaptive=True, density_window=50_000, adaptive_max_size=500_000
        )
        result = compute_intervals(variants, config)
        assert len(result) > 1
        for interval, group in result:
            assert interval.end - interval.start <= 500_000

    def test_interval_centered_on_variant(self):
        from hvantk.data.alphagenome_streamer import compute_intervals

        variants = [SimpleVariant("chr5", 1_000_000, "A", "T")]
        config = self._make_config(adaptive=False, default_size=200_000)
        result = compute_intervals(variants, config)
        interval = result[0][0]
        midpoint = (interval.start + interval.end) // 2
        assert abs(midpoint - 1_000_000) <= 1

    def test_empty_variant_list(self):
        from hvantk.data.alphagenome_streamer import compute_intervals

        result = compute_intervals([], self._make_config())
        assert result == []


class TestCheckpointManager:
    def test_fresh_start_no_checkpoints(self, tmp_path):
        from hvantk.data.alphagenome_streamer import CheckpointManager

        mgr = CheckpointManager(str(tmp_path / "out"))
        assert mgr.completed_intervals == []
        assert mgr.failed_variants == []

    def test_save_and_reload_state(self, tmp_path):
        from hvantk.data.alphagenome_streamer import CheckpointManager

        out_dir = str(tmp_path / "out")
        mgr = CheckpointManager(out_dir)
        mgr.mark_interval_complete("chr1:100000-200000")
        mgr.record_failed_variant("chr1", 150000, "A", "T", "timeout")
        mgr.save_state()

        mgr2 = CheckpointManager(out_dir)
        assert "chr1:100000-200000" in mgr2.completed_intervals
        assert len(mgr2.failed_variants) == 1
        assert mgr2.failed_variants[0]["chrom"] == "chr1"

    def test_save_batch_data(self, tmp_path):
        from hvantk.data.alphagenome_streamer import CheckpointManager

        out_dir = str(tmp_path / "out")
        mgr = CheckpointManager(out_dir)
        batch_data = {"variant_1": {"rna_seq": [1.0, 2.0]}}
        mgr.save_batch(0, batch_data)

        batch_path = os.path.join(out_dir, "_checkpoints", "batch_000.json")
        assert os.path.isfile(batch_path)
        with open(batch_path) as f:
            loaded = json.load(f)
        assert loaded == batch_data

    def test_is_interval_complete(self, tmp_path):
        from hvantk.data.alphagenome_streamer import CheckpointManager

        mgr = CheckpointManager(str(tmp_path / "out"))
        mgr.mark_interval_complete("chr1:100-200")
        assert mgr.is_interval_complete("chr1:100-200") is True
        assert mgr.is_interval_complete("chr2:100-200") is False

    def test_clear_checkpoints(self, tmp_path):
        from hvantk.data.alphagenome_streamer import CheckpointManager

        out_dir = str(tmp_path / "out")
        mgr = CheckpointManager(out_dir)
        mgr.mark_interval_complete("chr1:100-200")
        mgr.save_batch(0, {"data": True})
        mgr.save_state()
        mgr.clear()
        assert mgr.completed_intervals == []
        assert not os.path.isfile(
            os.path.join(out_dir, "_checkpoints", "state.json")
        )


class TestRateLimitedCaller:
    def _make_api_config(self, max_retries=3, retry_backoff=0.01, request_timeout=5):
        return {
            "api": {
                "key": "test-key",
                "max_retries": max_retries,
                "retry_backoff": retry_backoff,
                "request_timeout": request_timeout,
            }
        }

    def test_successful_call(self):
        from hvantk.data.alphagenome_streamer import RateLimitedCaller

        mock_model = MagicMock()
        mock_model.predict_variant.return_value = {"rna_seq": [1.0]}
        caller = RateLimitedCaller(mock_model, self._make_api_config())

        result = caller.call_predict_variant(
            interval=MagicMock(), variant=MagicMock(),
            ontology_terms=["UBERON:0001157"], output_types=["RNA_SEQ"],
        )
        assert result == {"rna_seq": [1.0]}
        mock_model.predict_variant.assert_called_once()

    def test_retry_on_transient_error(self):
        from hvantk.data.alphagenome_streamer import RateLimitedCaller

        mock_model = MagicMock()
        mock_model.predict_variant.side_effect = [
            Exception("503 Server Error"),
            Exception("503 Server Error"),
            {"rna_seq": [1.0]},
        ]
        caller = RateLimitedCaller(mock_model, self._make_api_config())
        result = caller.call_predict_variant(
            interval=MagicMock(), variant=MagicMock(),
            ontology_terms=[], output_types=[],
        )
        assert result == {"rna_seq": [1.0]}
        assert mock_model.predict_variant.call_count == 3

    def test_max_retries_exceeded_returns_none(self):
        from hvantk.data.alphagenome_streamer import RateLimitedCaller

        mock_model = MagicMock()
        mock_model.predict_variant.side_effect = Exception("503 Server Error")
        caller = RateLimitedCaller(
            mock_model, self._make_api_config(max_retries=2)
        )
        result = caller.call_predict_variant(
            interval=MagicMock(), variant=MagicMock(),
            ontology_terms=[], output_types=[],
        )
        assert result is None
        assert mock_model.predict_variant.call_count == 2

    def test_cooldown_after_consecutive_rate_limits(self):
        from hvantk.data.alphagenome_streamer import RateLimitedCaller

        mock_model = MagicMock()
        caller = RateLimitedCaller(mock_model, self._make_api_config())
        caller._consecutive_rate_limits = 3
        assert caller._should_cooldown() is True
        caller._consecutive_rate_limits = 2
        assert caller._should_cooldown() is False
