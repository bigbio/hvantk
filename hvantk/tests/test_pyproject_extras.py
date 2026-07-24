"""Guard the declared Poetry extras/deps that runtime imports rely on.

#198: cptac imports pyranges, whose compiled dep `sorted_nearest` is not always
pulled by an extras install, breaking `import cptac`. The ptm extra must declare
`sorted-nearest` explicitly.
"""
from __future__ import annotations

from pathlib import Path


def _load_pyproject() -> dict:
    try:
        import tomllib as toml  # Python >= 3.11
    except ModuleNotFoundError:  # Python 3.10
        import tomli as toml
    root = Path(__file__).resolve().parents[2]
    return toml.loads((root / "pyproject.toml").read_text())


def test_ptm_extra_includes_sorted_nearest():
    poetry = _load_pyproject()["tool"]["poetry"]
    assert "sorted-nearest" in poetry["extras"]["ptm"]
    assert "sorted-nearest" in poetry["dependencies"]
