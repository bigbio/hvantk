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


def test_lmm_module_imports_without_statsmodels_installed(monkeypatch):
    """#374 review item 4: `require_statsmodels()` used to run at MODULE scope
    (`smf = require_statsmodels()`), so merely `from hvantk.algorithms.ptm.lmm import
    run_lmm` failed on an install lacking statsmodels -- not only a call into the
    module. `require_scanpy`, the stated model for this pattern, is called inside the
    function that needs it; `lmm.py` must do the same.
    """
    _hide_statsmodels(monkeypatch)
    import hvantk.algorithms.ptm.lmm as lmm  # must not raise

    assert hasattr(lmm, "run_lmm")
    assert hasattr(lmm, "run_binned_interaction_lmm")


def test_run_lmm_raises_the_actionable_import_error_when_called_without_statsmodels(
    monkeypatch,
):
    import pandas as pd

    _hide_statsmodels(monkeypatch)
    import hvantk.algorithms.ptm.lmm as lmm

    df = pd.DataFrame(
        {"gene": ["G1"], "af_filled": [0.1], "is_ptm": [True]}
    )
    with pytest.raises(ImportError) as info:
        lmm.run_lmm(df, stratum="A")
    msg = str(info.value)
    assert "constraint" in msg and "hvantk[constraint]" in msg


def test_ptm_test_names_the_extra_instead_of_a_traceback(tmp_path, monkeypatch):
    from hvantk.tools.ptm.ptm_cli import ptm_group

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
