"""#362: a missing optional dependency must name its extra, not print a traceback."""
from __future__ import annotations

import sys

import pytest
from click.testing import CliRunner


def _hide_statsmodels(monkeypatch):
    for name in [m for m in list(sys.modules) if m == "statsmodels" or m.startswith("statsmodels.")]:
        monkeypatch.delitem(sys.modules, name)
    monkeypatch.setitem(sys.modules, "statsmodels", None)  # `import statsmodels` -> ModuleNotFoundError
    monkeypatch.delitem(sys.modules, "hvantk.algorithms.ptm.lmm", raising=False)


def test_require_statsmodels_names_the_constraint_extra(monkeypatch):
    from hvantk.algorithms.ptm.optional_deps import require_statsmodels

    _hide_statsmodels(monkeypatch)
    with pytest.raises(ImportError) as info:
        require_statsmodels()
    msg = str(info.value)
    assert "constraint" in msg and "hvantk[constraint]" in msg and "poetry install --extras constraint" in msg


def test_require_statsmodels_returns_the_formula_api_when_installed():
    pytest.importorskip("statsmodels")
    from hvantk.algorithms.ptm.optional_deps import require_statsmodels

    assert hasattr(require_statsmodels(), "mixedlm")


def test_ptm_test_names_the_extra_instead_of_a_traceback(tmp_path, monkeypatch):
    from hvantk.tools.ptm.ptm_cli import ptm_group  # adjust to the group/command object the existing ptm CLI tests use

    _hide_statsmodels(monkeypatch)
    inp = tmp_path / "in.tsv"
    inp.write_text("gene\taf\tis_ptm\tstratum\nG1\t0.1\tTrue\tA\nG2\t0.2\tFalse\tA\n")
    result = CliRunner().invoke(
        ptm_group,
        ["test", "--test", "lmm", "--input", str(inp), "--stratum-col", "stratum",
         "--output", str(tmp_path / "out.tsv")],
    )
    assert result.exit_code != 0
    assert "constraint" in result.output, result.output
    assert "Traceback" not in result.output
    assert "No module named" not in result.output
