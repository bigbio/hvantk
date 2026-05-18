import textwrap
from pathlib import Path
from click.testing import CliRunner
import yaml

from hvantk.tools.infra.catalog_cli import catalog


def _write_catalog(tmpdir: Path, content: str) -> Path:
    p = tmpdir / "catalog.yaml"
    p.write_text(textwrap.dedent(content))
    return p


def test_catalog_list_on_packaged(monkeypatch, tmp_path):
    # Create list style catalog
    cat_path = _write_catalog(
        tmp_path,
        """
    - id: test.entry
      provenance:
        builder: hvantk.tests.test_cli_catalog._dummy_builder
      params:
        answer: 42
    """,
    )

    runner = CliRunner()
    result = runner.invoke(catalog, ["list", "--catalog", str(cat_path)])
    assert result.exit_code == 0
    assert "test.entry" in result.output


def _dummy_builder(answer=0):  # pragma: no cover - executed via CLI
    return f"answer={answer}"


def test_catalog_show(monkeypatch, tmp_path):
    cat_path = _write_catalog(
        tmp_path,
        """
    - id: demo
      provenance:
        builder: hvantk.tests.test_cli_catalog._dummy_builder
      params:
        foo: bar
    """,
    )
    runner = CliRunner()
    result = runner.invoke(catalog, ["show", "demo", "--catalog", str(cat_path)])
    assert result.exit_code == 0
    assert "foo: bar" in result.output


def test_catalog_build(monkeypatch, tmp_path):
    cat_path = _write_catalog(
        tmp_path,
        """
    - id: buildme
      provenance:
        builder: hvantk.tests.test_cli_catalog._dummy_builder
      params:
        answer: 5
    """,
    )
    runner = CliRunner()
    result = runner.invoke(
        catalog,
        ["build", "buildme", "--catalog", str(cat_path), "--override", "answer=7"],
    )
    assert result.exit_code == 0
    assert "Invoking builder" in result.output
    assert "answer=7" in result.output
    assert "Built buildme" in result.output
