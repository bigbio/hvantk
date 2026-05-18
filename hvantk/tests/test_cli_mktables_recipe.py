import json
from click.testing import CliRunner
from unittest.mock import patch

from hvantk.tools.build.make_table_batch_cli import mktable_batch_cli


def test_mktables_recipe_success(tmp_path):
    recipe = {
        "tables": [
            {
                "name": "clinvar:variants",
                "input": "/data/clin.vcf.bgz",
                "output": "/out/clinvar.ht",
                "params": {"reference_genome": "GRCh38", "export_tsv": True},
            },
            {
                "name": "gevir",
                "input": "/data/gevir.tsv.bgz",
                "output": "/out/gevir.ht",
                "params": {"fields": "oe_syn_upper,oe_mis_upper"},
            },
        ]
    }
    recipe_path = tmp_path / "tables.json"
    recipe_path.write_text(json.dumps(recipe))

    runner = CliRunner()
    with patch("hvantk.tools.build.make_table_batch_cli._load_recipe") as mock_load, patch(
        "hvantk.tools.build.make_table_batch_cli.logger"
    ), patch("hvantk.tables.registry.run_table_builder") as mock_runner:
        mock_load.return_value = recipe
        result = runner.invoke(mktable_batch_cli, ["--recipe", str(recipe_path)])
        assert result.exit_code == 0
        # Two calls
        assert mock_runner.call_count == 2
        mock_runner.assert_any_call(
            "clinvar:variants",
            "/data/clin.vcf.bgz",
            "/out/clinvar.ht",
            {"reference_genome": "GRCh38", "export_tsv": True},
        )
        mock_runner.assert_any_call(
            "gevir",
            "/data/gevir.tsv.bgz",
            "/out/gevir.ht",
            {"fields": "oe_syn_upper,oe_mis_upper"},
        )
        assert "table created" in result.output


def test_mktables_recipe_missing_tables(tmp_path):
    recipe = {}
    recipe_path = tmp_path / "empty.json"
    recipe_path.write_text(json.dumps(recipe))

    runner = CliRunner()
    result = runner.invoke(mktable_batch_cli, ["--recipe", str(recipe_path)])
    assert result.exit_code != 0
    assert "Recipe is missing non-empty 'tables' list" in result.output


def test_mktables_recipe_bad_entry(tmp_path):
    recipe = {"tables": [{"name": "clinvar", "input": "/x"}]}  # missing output
    recipe_path = tmp_path / "bad.json"
    recipe_path.write_text(json.dumps(recipe))

    runner = CliRunner()
    result = runner.invoke(mktable_batch_cli, ["--recipe", str(recipe_path)])
    assert result.exit_code != 0
    assert "Invalid entry in recipe" in result.output
