# hvantk.core package
"""Shared infrastructure layer.

The four ``backends`` symbols below stay importable from ``hvantk.core`` for
backwards compatibility, but are resolved lazily (PEP 562). Importing them
eagerly pulled in ``hvantk.core.models.__init__`` -- and with it pandas and the
concrete artifact classes -- on *every* ``hvantk.core.*`` import, which cost
~1.7 s before any work was done. Nothing in the tree actually imports these
names from here; they are kept only so external callers do not break.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - import-time typing only
    from hvantk.core.models.backends import (
        AlgorithmMeta,
        Backend,
        algorithm,
        get_algorithm_meta,
    )

__all__ = [
    "AlgorithmMeta",
    "Backend",
    "algorithm",
    "get_algorithm_meta",
]


def __getattr__(name: str):
    if name in __all__:
        import importlib

        value = getattr(importlib.import_module("hvantk.core.models.backends"), name)
        globals()[name] = value  # cache, so the import cost is paid once
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    # Union, not just __all__: returning only the lazy names would hide the
    # real subpackages (config, io, models, plugin, utils, ...) from dir(),
    # REPL completion and inspect.getmembers even after they are imported.
    return sorted({*__all__, *globals()} - {"TYPE_CHECKING"})
