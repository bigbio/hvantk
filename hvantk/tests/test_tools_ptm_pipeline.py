"""Tests for the PTM workflow wrapper in tools/ptm/."""

from __future__ import annotations

import pytest


def test_tools_ptm_exports_download_and_pipeline():
    from hvantk.tools.ptm.pipeline import (
        download_uniprot_ptm,
        ptm_build_pipeline,
    )

    assert callable(download_uniprot_ptm)
    assert callable(ptm_build_pipeline)


def test_algorithms_ptm_has_no_skill_imports():
    """AST-walk hvantk/algorithms/ptm/ asserting no hvantk.skills import.

    This is a tighter check than test_dependency_directions; it specifically
    catches any regression in the ptm subpackage.
    """
    import ast
    from pathlib import Path

    root = Path(__file__).resolve().parents[1] / "algorithms" / "ptm"
    bad = []
    for py in root.rglob("*.py"):
        if "__pycache__" in py.parts:
            continue
        tree = ast.parse(py.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name.startswith("hvantk.skills"):
                        bad.append((py.name, alias.name))
            elif isinstance(node, ast.ImportFrom) and node.module:
                if node.module.startswith("hvantk.skills"):
                    bad.append((py.name, node.module))
    assert not bad, f"algorithms/ptm/ has skill imports: {bad}"


def test_ptm_build_pipeline_resolves_hyphenated_uniprot_key(monkeypatch, tmp_path):
    """#198 regression: the Hail build step must resolve 'uniprot-ptm:sites'
    (hyphen, per plugin.yaml `name: uniprot-ptm`), not the underscore form.

    The bug wasted the whole expensive GTF download + coordinate mapping run
    before KeyError-ing at the final build step.
    """
    import hvantk.tools.ptm.pipeline as pl
    from hvantk.algorithms.ptm.pipeline import PTMBuildConfig

    captured = {}

    class _FakeSpec:
        plugin_version = "0.1.0"

    class _FakeRegistry:
        def get_dataset(self, key):
            captured["key"] = key
            return _FakeSpec()

    class _FakeResult:
        n_mapped = 5
        # The combined name, so a builder fed the UniProt-only fallback path fails.
        mapped_tsv_path = str(tmp_path / "ptm_sites_combined.tsv.bgz")
        output_ht = None

    # get_registry + run_builder_for_spec are imported INSIDE ptm_build_pipeline;
    # patch them at their source modules. ptm_build_pipeline_core is a module global.
    monkeypatch.setattr(
        "hvantk.core.plugin.loader.get_registry", lambda: _FakeRegistry()
    )
    monkeypatch.setattr(
        "hvantk.core.plugin.run_builder.run_builder_for_spec",
        lambda *a, **k: captured.update(k),
    )
    monkeypatch.setattr(pl, "ptm_build_pipeline_core", lambda cfg: _FakeResult())

    # The file must exist: ptm_build_pipeline validates the config before it
    # downloads or maps anything.
    ptm_tsv = tmp_path / "ptm.tsv"
    ptm_tsv.write_text("")
    cfg = PTMBuildConfig(
        output_dir=str(tmp_path),
        output_ht=str(tmp_path / "out.ht"),
        ptm_tsv=str(ptm_tsv),  # non-None -> download step skipped
    )
    result = pl.ptm_build_pipeline(cfg)

    assert captured["key"] == "uniprot-ptm:sites"
    # A multi-source run must build the table from the combined TSV, and say so.
    assert captured["parsed_input"] == _FakeResult.mapped_tsv_path
    assert result.output_ht == cfg.output_ht


def test_ptm_build_pipeline_fails_when_no_site_maps(monkeypatch, tmp_path):
    """A build that maps no PTM site must raise before the Hail build step.

    Returning without a table would let the CLI exit 0 while a table from an
    earlier run at ``output_ht`` stays there for ``ptm annotate`` to read.
    """
    import hvantk.tools.ptm.pipeline as pl
    from hvantk.algorithms.ptm.pipeline import PTMBuildConfig, PTMBuildResult

    class _FakeSpec:
        plugin_version = "0.1.0"

    class _FakeRegistry:
        def get_dataset(self, key):
            return _FakeSpec()

    builds = []
    monkeypatch.setattr(
        "hvantk.core.plugin.loader.get_registry", lambda: _FakeRegistry()
    )
    monkeypatch.setattr(
        "hvantk.core.plugin.run_builder.run_builder_for_spec",
        lambda *a, **k: builds.append(k),
    )
    monkeypatch.setattr(
        pl, "ptm_build_pipeline_core", lambda cfg: PTMBuildResult(n_total=3, n_mapped=0)
    )

    ptm_tsv = tmp_path / "ptm.tsv"
    ptm_tsv.write_text("")
    earlier_table = tmp_path / "out.ht"
    earlier_table.mkdir()
    (earlier_table / "marker").write_text("earlier run")
    cfg = PTMBuildConfig(
        output_dir=str(tmp_path),
        output_ht=str(earlier_table),
        ptm_tsv=str(ptm_tsv),
    )

    with pytest.raises(ValueError, match="No PTM sites mapped"):
        pl.ptm_build_pipeline(cfg)
    assert builds == []
    assert (earlier_table / "marker").read_text() == "earlier run"
