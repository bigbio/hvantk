"""The PEP 562 installer shared by the package inits that must stay cheap to import (#306)."""
from __future__ import annotations

import sys
import types

import pytest

from hvantk.core.utils.lazy_exports import install_lazy_exports


def _package(name: str, exports: dict, **kwargs) -> types.ModuleType:
    pkg = types.ModuleType(name)
    pkg.__path__ = []  # marks it as a package
    sys.modules[name] = pkg
    install_lazy_exports(pkg.__dict__, exports, **kwargs)
    return pkg


def test_export_is_resolved_on_first_access_and_cached(monkeypatch):
    pkg = _package("lazypkg_a", {"sqrt": "math"})
    monkeypatch.delitem(sys.modules, "lazypkg_a", raising=False)
    assert "sqrt" not in pkg.__dict__
    import math

    assert pkg.sqrt is math.sqrt
    assert pkg.__dict__["sqrt"] is math.sqrt  # cached, so the import cost is paid once
    assert "sqrt" in dir(pkg)


def test_missing_dependency_gets_the_hint(monkeypatch):
    pkg = _package(
        "lazypkg_b",
        {"thing": "definitely_not_installed_xyz"},
        missing_hint=lambda exc: f"install the extra for {exc.name}",
    )
    monkeypatch.delitem(sys.modules, "lazypkg_b", raising=False)
    with pytest.raises(ImportError, match="install the extra for definitely_not_installed_xyz"):
        pkg.thing


def test_unknown_attribute_raises_attribute_error(monkeypatch):
    pkg = _package("lazypkg_c", {})
    monkeypatch.delitem(sys.modules, "lazypkg_c", raising=False)
    with pytest.raises(AttributeError, match="lazypkg_c.*nope"):
        pkg.nope


def test_submodule_access_still_works_for_a_real_package():
    """`import hvantk.algorithms.enrichex as ex; ex.constants` must not become an
    AttributeError once the package is lazy."""
    import hvantk.algorithms.enrichex as ex

    assert ex.constants.GENOTYPE_AGGREGATION_METHODS
