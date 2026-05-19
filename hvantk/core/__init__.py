# hvantk.core package
from hvantk.core.models.backends import (
    AlgorithmMeta,
    Backend,
    algorithm,
    get_algorithm_meta,
)
from hvantk.core.utils.readers import (
    DataReader,
    DuckDBReader,
    HailReader,
    PandasReader,
)
from hvantk.core.utils.router import BackendRouter, ReaderFactory
from hvantk.core.utils.writers import HailTableWriter

__all__ = [
    "AlgorithmMeta",
    "Backend",
    "BackendRouter",
    "DataReader",
    "DuckDBReader",
    "HailReader",
    "HailTableWriter",
    "PandasReader",
    "ReaderFactory",
    "algorithm",
    "get_algorithm_meta",
]

# ---------------------------------------------------------------------------
# Backwards-compat shims for callers using the pre-restructure import paths.
# REMOVE IN PHASE 7 of the restructure plan.
# ---------------------------------------------------------------------------
import warnings as _warnings


def __getattr__(name: str):
    """Forward old `from hvantk.core import plugin_loader` style imports to
    their new sub-package locations, with a DeprecationWarning so callers
    can find the import to update.
    """
    _shim_map = {
        # plugin runtime
        "plugin_api":              "hvantk.core.plugin.api",
        "plugin_loader":           "hvantk.core.plugin.loader",
        "drift_runner":            "hvantk.core.plugin.drift_runner",
        # tool runtime
        "tool_api":                "hvantk.core.tool.api",
        "tool_loader":             "hvantk.core.tool.loader",
        # models
        "anndata_utils":           "hvantk.core.models.anndata_utils",
        "backends":                "hvantk.core.models.backends",
        "metadata":                "hvantk.core.models.metadata",
        # utils
        "bgzf":                    "hvantk.core.utils.bgzf",
        "converters":              "hvantk.core.utils.converters",
        "readers":                 "hvantk.core.utils.readers",
        "writers":                 "hvantk.core.utils.writers",
        "router":                  "hvantk.core.utils.router",
        "hail_context":            "hvantk.core.utils.hail_context",
    }
    if name in _shim_map:
        _warnings.warn(
            f"hvantk.core.{name} moved to {_shim_map[name]}; "
            "update the import (this shim is removed in restructure Phase 7)",
            DeprecationWarning,
            stacklevel=2,
        )
        import importlib
        return importlib.import_module(_shim_map[name])
    raise AttributeError(f"module 'hvantk.core' has no attribute '{name}'")
