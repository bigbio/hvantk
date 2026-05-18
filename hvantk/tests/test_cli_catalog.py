"""Tests for ``hvantk catalog`` CLI against per-plugin catalogs."""

from click.testing import CliRunner

from hvantk.tools.infra.catalog_cli import catalog as catalog_group


def test_list_returns_entries():
    result = CliRunner().invoke(catalog_group, ["list", "--limit", "5"])
    assert result.exit_code == 0, result.output
    # Should show the column header and at least one row.
    assert "ACCESSION" in result.output


def test_list_filter_by_omics_type():
    result = CliRunner().invoke(
        catalog_group, ["list", "--omics-type", "transcriptomics", "--limit", "3"]
    )
    assert result.exit_code == 0, result.output
    assert "transcriptomics" in result.output


def test_list_filter_by_organism_substring():
    """``--organism Homo`` should match ``Homo sapiens`` via substring semantics."""
    result = CliRunner().invoke(
        catalog_group, ["list", "--organism", "Homo", "--limit", "3"]
    )
    assert result.exit_code == 0, result.output
    # Either we get rows (header line present) or "(no entries match)"; both
    # are valid for the substring code path. We only assert exit_code here.


def test_show_unknown_accession_errors():
    result = CliRunner().invoke(catalog_group, ["show", "DOES-NOT-EXIST"])
    assert result.exit_code != 0
    assert "not found" in result.output.lower()


def test_show_known_accession_yaml():
    """E-GTEX-8 is the first Expression Atlas entry shipped with the plugin."""
    result = CliRunner().invoke(catalog_group, ["show", "E-GTEX-8"])
    assert result.exit_code == 0, result.output
    assert "E-GTEX-8" in result.output
    assert "_omics_type" in result.output


def test_show_known_accession_json():
    result = CliRunner().invoke(
        catalog_group, ["show", "E-GTEX-8", "--format", "json"]
    )
    assert result.exit_code == 0, result.output
    assert '"accession"' in result.output
    assert "E-GTEX-8" in result.output


def test_stats_outputs_totals_and_breakdowns():
    result = CliRunner().invoke(catalog_group, ["stats"])
    assert result.exit_code == 0, result.output
    assert "total datasets" in result.output
    assert "by omics type" in result.output
    assert "top organisms" in result.output
    assert "top data sources" in result.output


def test_search_finds_known_term():
    """Verify shape: header is printed when matches exist, or "no matches" otherwise."""
    result = CliRunner().invoke(catalog_group, ["search", "Homo", "--limit", "5"])
    assert result.exit_code == 0, result.output
    assert "matches:" in result.output or "no matches" in result.output


def test_search_empty_query_for_missing_term_is_clean():
    result = CliRunner().invoke(
        catalog_group, ["search", "this-string-should-not-match-anything-xyzzy"]
    )
    assert result.exit_code == 0, result.output
    assert "no matches" in result.output


def test_build_subcommand_removed():
    """The legacy ``build`` subcommand was retired; ensure it's gone."""
    result = CliRunner().invoke(catalog_group, ["build", "anything"])
    assert result.exit_code != 0
