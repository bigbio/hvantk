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


@pytest.mark.parametrize(
    "package_name",
    [
        "hvantk.algorithms.enrichex",
        "hvantk.algorithms.qtlcascade",
    ],
)
def test_all_lazy_exports_are_resolvable(package_name):
    """Every name in a lazy package's __all__ must resolve on access.

    Wrong module paths in the exports map are now only caught at attribute access,
    not at package import time. This test ensures a typo in the map does not silently
    break the export: __dir__ lists every map key, so set(__all__) - set(dir())
    cannot catch it. We must actually getattr() each name.
    """
    pkg = __import__(package_name, fromlist=["__all__"])
    exported = set(pkg.__all__)
    unresolvable = []
    for name in sorted(exported):
        try:
            getattr(pkg, name)
        except Exception as exc:  # noqa: BLE001 - report any failure
            unresolvable.append(f"{name}: {type(exc).__name__}: {exc}")
    assert (
        not unresolvable
    ), f"{package_name} lists these in __all__ but cannot resolve them: {unresolvable}"
