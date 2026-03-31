import os
from dataclasses import dataclass

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
