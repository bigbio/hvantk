"""Smoke test: the ptm workflow wrapper exists in tools/ptm/ after Phase C."""
from __future__ import annotations


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
        mapped_tsv_path = str(tmp_path / "ptm_sites_mapped.tsv.bgz")
        output_ht = None

    # get_registry + run_builder_for_spec are imported INSIDE ptm_build_pipeline;
    # patch them at their source modules. ptm_build_pipeline_core is a module global.
    monkeypatch.setattr(
        "hvantk.core.plugin.loader.get_registry", lambda: _FakeRegistry()
    )
    monkeypatch.setattr(
        "hvantk.core.plugin.run_builder.run_builder_for_spec", lambda *a, **k: None
    )
    monkeypatch.setattr(pl, "ptm_build_pipeline_core", lambda cfg: _FakeResult())

    cfg = PTMBuildConfig(
        output_dir=str(tmp_path),
        output_ht=str(tmp_path / "out.ht"),
        ptm_tsv=str(tmp_path / "ptm.tsv"),  # non-None -> download step skipped
    )
    pl.ptm_build_pipeline(cfg)

    assert captured["key"] == "uniprot-ptm:sites"
