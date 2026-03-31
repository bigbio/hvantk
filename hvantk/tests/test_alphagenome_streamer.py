import os
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
