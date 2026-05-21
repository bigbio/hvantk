"""
Registry and adapters for table builders.

Allows running named builders with a unified interface from recipes.

Contract per entry:
- name: str (e.g., "clinvar", "insider:variants", "gevir", "gnomad-metrics", "ensembl-gene")
- input_path: str
- output_path: str
- params: dict (optional) – builder-specific parameters

This module adapts hvantk.core.builders.table functions to this contract.

Internal dispatch tables
------------------------
``_TABLE_BUILDERS`` and ``_MATRIX_BUILDERS`` are module-private implementation
details.  External callers should use the public dispatch functions
``run_table_builder`` / ``run_matrix_builder`` instead of importing the dicts
directly.
"""

from __future__ import annotations

import logging
import inspect
from typing import Callable, Dict, Any, get_type_hints

logger = logging.getLogger(__name__)


# Helper functions for parameter type conversion


def _parse_list_param(value: Any) -> list | None:
    """Parse a generic list parameter (comma-separated string or list)."""
    if value is None:
        return None
    if isinstance(value, list):
        return value
    if isinstance(value, str):
        parts = [x.strip() for x in value.split(",") if x.strip()]
        return parts or None
    return None


def _convert_param_type(param_name: str, value: Any, target_type: type) -> Any:
    """
    Convert a parameter value to the expected type based on function signature.

    Handles common conversions:
    - bool: converts strings and integers
    - int: parses numeric values
    - str: ensures string type
    - list: handles comma-separated strings (for fields, prefixes, etc.)
    """
    if value is None:
        return None

    # Special handling for list-like parameters
    if param_name in ("fields", "group_prefixes", "categorical_cols", "numeric_cols"):
        return _parse_list_param(value)

    # Type-specific conversions
    if target_type is bool:
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            return value.lower() in ("true", "1", "yes", "on")
        return bool(value)
    elif target_type is int:
        return int(value)
    elif target_type is str:
        return str(value)
    elif target_type is list:
        return _parse_list_param(value)

    return value


def create_table_adapter(
    builder_module: str, builder_function: str
) -> Callable[[str, str, Dict[str, Any] | None], None]:
    """
    Factory function to create adapter functions for table builders.

    Uses introspection to automatically map recipe params to builder function parameters.
    This eliminates the need for manual adapter functions.

    Parameters
    ----------
    builder_module : str
        The module path containing the builder function (e.g., 'hvantk.core.builders.table').
    builder_function : str
        The name of the builder function (e.g., 'create_clinvar_tb').

    Returns
    -------
    Callable
        An adapter function with signature: (input_path, output_path, params) -> None

    Example
    -------
    >>> clinvar_adapter = create_table_adapter(
    ...     'hvantk.core.builders.table',
    ...     'create_clinvar_tb'
    ... )

    """

    def adapter(
        input_path: str, output_path: str, params: Dict[str, Any] | None = None
    ):
        # Lazy import to avoid importing Hail at CLI import time
        import importlib

        module = importlib.import_module(builder_module)
        builder_func = getattr(module, builder_function)

        # Get function signature for parameter mapping
        sig = inspect.signature(builder_func)
        params = params or {}

        # Get resolved type hints (handles string annotations from __future__)
        try:
            type_hints = get_type_hints(builder_func)
        except Exception:
            # If get_type_hints fails, fall back to raw annotations
            type_hints = {}

        # Build kwargs for the builder function
        kwargs = {"input_path": input_path, "output_path": output_path}

        # Map params to function parameters with type conversion
        for param_name, param_obj in sig.parameters.items():
            if param_name in ("input_path", "output_path"):
                continue  # Already handled

            if param_name in params:
                value = params[param_name]
                # Get target type from type hints, annotation, or default value
                target_type = None
                if param_name in type_hints:
                    target_type = type_hints[param_name]
                elif param_obj.annotation != inspect.Parameter.empty:
                    target_type = param_obj.annotation
                elif param_obj.default != inspect.Parameter.empty:
                    target_type = type(param_obj.default)

                # Convert parameter value to expected type
                if target_type:
                    value = _convert_param_type(param_name, value, target_type)

                kwargs[param_name] = value

        # Call the builder function
        builder_func(**kwargs)

    return adapter


def create_matrix_adapter(
    builder_module: str, builder_function: str, required_inputs: list[str] | None = None
) -> Callable[[Dict[str, str], str, Dict[str, Any] | None], None]:
    """
    Factory function to create adapter functions for matrix builders.

    Similar to create_table_adapter but handles multi-input matrix builders.

    Parameters
    ----------
    builder_module : str
        The module path containing the builder function.
    builder_function : str
        The name of the builder function.
    required_inputs : list of str, optional
        List of required input keys that must be present in the inputs dict.

    Returns
    -------
    Callable
        An adapter function with signature: (inputs, output_mt, params) -> None

    """

    def adapter(
        inputs: Dict[str, str], output_mt: str, params: Dict[str, Any] | None = None
    ):
        # Validate required inputs
        if required_inputs:
            missing = [k for k in required_inputs if k not in inputs or not inputs[k]]
            if missing:
                raise ValueError(
                    f"Missing required inputs for '{builder_function}': {missing}"
                )

        # Lazy import
        import importlib

        module = importlib.import_module(builder_module)
        builder_func = getattr(module, builder_function)

        # Get function signature
        sig = inspect.signature(builder_func)
        params = params or {}

        # Get resolved type hints (handles string annotations from __future__)
        try:
            type_hints = get_type_hints(builder_func)
        except Exception:
            # If get_type_hints fails, fall back to raw annotations
            type_hints = {}

        # Build kwargs - map inputs to function parameters
        kwargs = {}

        # Map input paths to function parameters
        for param_name, param_obj in sig.parameters.items():
            # Check if this parameter corresponds to an input file
            if param_name.endswith("_path") or param_name.endswith("_file"):
                # Try to find matching input key
                for input_key, input_path in inputs.items():
                    # Match based on naming conventions
                    if (
                        input_key in param_name
                        or param_name.replace("_path", "").replace("_file", "")
                        in input_key
                    ):
                        kwargs[param_name] = input_path
                        break
            elif param_name in ("output_mt", "output_path"):
                kwargs[param_name] = output_mt
            elif param_name in params:
                # Handle parameters from params dict
                value = params[param_name]
                # Get target type from type hints, annotation, or default value
                target_type = None
                if param_name in type_hints:
                    target_type = type_hints[param_name]
                elif param_obj.annotation != inspect.Parameter.empty:
                    target_type = param_obj.annotation
                elif param_obj.default != inspect.Parameter.empty:
                    target_type = type(param_obj.default)

                if target_type:
                    value = _convert_param_type(param_name, value, target_type)

                kwargs[param_name] = value

        builder_func(**kwargs)

    return adapter


# Table builders — populated by the plugin loader via _initialize_plugin_registrations.
# All 7 hardcoded entries (gevir, gnomad-metrics, ensembl-gene, dbnsfp, cosmic-cgc,
# pqtl, alphagenome) were promoted to Phase B plugins in Phase K and removed here.
# Backward-compat short-name aliases are added by _add_phase_k_aliases() after the
# plugin loader runs; see that function below.
_TABLE_BUILDERS: Dict[str, Callable[[str, str, Dict[str, Any] | None], None]] = {}


def run_table_builder(
    name: str, input_path: str, output_path: str, params: Dict[str, Any] | None = None
) -> None:
    """Run a registered table builder by name.

    Raises KeyError if the builder name is unknown.
    """
    if name not in _TABLE_BUILDERS:
        raise KeyError(f"Unknown table builder: {name}")
    logger.info(
        f"Running builder '{name}' with input={input_path} output={output_path} params={params}"
    )
    _TABLE_BUILDERS[name](input_path, output_path, params or {})


# Matrix builders - automatically generated adapters using factory
# All entries now come from the plugin registry via
# ``_initialize_plugin_registrations`` below.
_MATRIX_BUILDERS: Dict[
    str, Callable[[Dict[str, str], str, Dict[str, Any] | None], None]
] = {}


def run_matrix_builder(
    name: str,
    inputs: Dict[str, str],
    output_mt: str,
    params: Dict[str, Any] | None = None,
) -> None:
    """Run a registered matrix builder by name.

    Raises KeyError if the builder name is unknown.
    """
    if name not in _MATRIX_BUILDERS:
        raise KeyError(f"Unknown matrix builder: {name}")
    logger.info(
        f"Running matrix builder '{name}' with inputs={inputs} output={output_mt} params={params}"
    )
    _MATRIX_BUILDERS[name](inputs, output_mt, params or {})


# --- Plugin-driven registrations (added by feat/data-handlers-refactoring) ---


def _apply_plugin_registrations(reg) -> None:
    """Add plugin-discovered builders to the internal dispatch dicts.

    Coexistence: this runs ALONGSIDE the create_table_adapter() block above.
    As each provider migrates to the plugin layout in follow-up plans, its
    create_table_adapter line is removed in the same migration commit; this
    function continues to populate the dict from the plugin.
    """
    from hvantk.core.plugin.api import DatasetSpec  # local import to avoid cycle

    def _wrap_builder(spec: "DatasetSpec"):
        if spec.backend == "hail":
            def adapter(input_path, output_path, params=None):
                spec.builder(input_path, output_path, **(params or {}))
            return adapter
        if spec.backend == "anndata":
            def adapter(inputs, output_mt, params=None):
                # AnnData builders take multi-input dicts — Phase 0 simplification:
                # the input dict is expanded directly as kwargs. New plugin builders
                # MUST name their kwargs to match the input dict keys.
                spec.builder(**inputs, output_path=output_mt, **(params or {}))
            return adapter
        # pandas backend: same signature as hail for now.
        def adapter(input_path, output_path, params=None):
            spec.builder(input_path, output_path, **(params or {}))
        return adapter

    for ds in reg.list_datasets(backend="hail"):
        _TABLE_BUILDERS[ds.name] = _wrap_builder(ds)
    for ds in reg.list_datasets(backend="pandas"):
        _TABLE_BUILDERS[ds.name] = _wrap_builder(ds)
    for ds in reg.list_datasets(backend="anndata"):
        _MATRIX_BUILDERS[ds.name] = _wrap_builder(ds)


def _add_phase_k_aliases() -> None:
    """Map the old short names to the new plugin-loaded compound keys.

    Preserves backward compatibility with existing mktable-batch recipes and
    CLI invocations that used the pre-Phase-K short names (gevir, gnomad-metrics,
    etc.). Runs once after _apply_plugin_registrations populates the long names.
    """
    aliases = {
        "gevir": "gevir:metrics",
        "gnomad-metrics": "gnomad-metrics:metrics",
        "ensembl-gene": "ensembl-gene:genes",
        "dbnsfp": "dbnsfp:variants",
        "cosmic-cgc": "cosmic-cgc:submissions",
        "pqtl": "pqtl:metrics",
        "alphagenome": "alphagenome:predictions",
    }
    for old, new in aliases.items():
        if new in _TABLE_BUILDERS and old not in _TABLE_BUILDERS:
            _TABLE_BUILDERS[old] = _TABLE_BUILDERS[new]


def _initialize_plugin_registrations() -> None:
    """Wire the module-level PluginRegistry into _TABLE_BUILDERS / _MATRIX_BUILDERS.

    Called once at module import time. Safe to call repeatedly; subsequent
    calls are no-ops because the registry is a module-level singleton.
    """
    from hvantk.core.plugin import loader as plugin_loader

    try:
        reg = plugin_loader.get_registry()
    except Exception as exc:  # noqa: BLE001 — never let plugin failure break hvantk import
        logger.warning("plugin loader failed to initialize: %s", exc)
        return
    _apply_plugin_registrations(reg)
    _add_phase_k_aliases()


_initialize_plugin_registrations()
