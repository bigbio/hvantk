"""
Registry and adapters for table builders.

Allows running named builders with a unified interface from recipes.

Contract per entry:
- name: str (e.g., "clinvar", "interactome", "gevir", "gnomad-metrics", "ensembl-gene")
- input_path: str
- output_path: str
- params: dict (optional) – builder-specific parameters

This module adapts hvantk.tables.table_builders functions to this contract.
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
        The module path containing the builder function (e.g., 'hvantk.tables.table_builders').
    builder_function : str
        The name of the builder function (e.g., 'create_clinvar_tb').

    Returns
    -------
    Callable
        An adapter function with signature: (input_path, output_path, params) -> None

    Example
    -------
    >>> clinvar_adapter = create_table_adapter(
    ...     'hvantk.tables.table_builders',
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


# Table builders - automatically generated adapters using factory
TABLE_BUILDERS: Dict[str, Callable[[str, str, Dict[str, Any] | None], None]] = {
    "clinvar": create_table_adapter(
        "hvantk.tables.table_builders", "create_clinvar_tb"
    ),
    "interactome": create_table_adapter(
        "hvantk.tables.table_builders", "create_interactome_tb"
    ),
    "gevir": create_table_adapter("hvantk.tables.table_builders", "create_gevir_tb"),
    "gnomad-metrics": create_table_adapter(
        "hvantk.tables.table_builders", "create_gnomad_constraint_gene_metrics_tb"
    ),
    "ensembl-gene": create_table_adapter(
        "hvantk.tables.table_builders", "create_ensembl_gene_tb"
    ),
    "dbnsfp": create_table_adapter("hvantk.tables.table_builders", "create_dbnsfp_tb"),
    "clingen-gene-disease": create_table_adapter(
        "hvantk.tables.table_builders", "create_clingen_gene_disease_tb"
    ),
    "gencc-submissions": create_table_adapter(
        "hvantk.tables.table_builders", "create_gencc_submissions_tb"
    ),
    "cosmic-cgc": create_table_adapter(
        "hvantk.tables.table_builders", "create_cosmic_cgc_tb"
    ),
    "ptm-sites": create_table_adapter(
        "hvantk.tables.table_builders", "create_ptm_sites_tb"
    ),
    "eqtl": create_table_adapter(
        "hvantk.tables.table_builders", "create_eqtl_tb"
    ),
    "pqtl": create_table_adapter(
        "hvantk.tables.table_builders", "create_pqtl_tb"
    ),
}


def run_table_builder(
    name: str, input_path: str, output_path: str, params: Dict[str, Any] | None = None
) -> None:
    """Run a registered table builder by name.

    Raises KeyError if the builder name is unknown.
    """
    if name not in TABLE_BUILDERS:
        raise KeyError(f"Unknown table builder: {name}")
    logger.info(
        f"Running builder '{name}' with input={input_path} output={output_path} params={params}"
    )
    TABLE_BUILDERS[name](input_path, output_path, params or {})


# Matrix builders - automatically generated adapters using factory
MATRIX_BUILDERS: Dict[
    str, Callable[[Dict[str, str], str, Dict[str, Any] | None], None]
] = {
    "ucsc": create_matrix_adapter(
        "hvantk.tables.matrix_builders",
        "build_ucsc_mt",
        required_inputs=["expression_matrix", "metadata"],
    ),
    "expression-atlas": create_matrix_adapter(
        "hvantk.tables.matrix_builders",
        "build_expression_atlas_mt",
        required_inputs=["expression_matrix", "sdrf"],
    ),
    "cptac": create_matrix_adapter(
        "hvantk.tables.matrix_builders",
        "build_cptac_mt",
        required_inputs=["expression", "metadata"],
    ),
}


def run_matrix_builder(
    name: str,
    inputs: Dict[str, str],
    output_mt: str,
    params: Dict[str, Any] | None = None,
) -> None:
    """Run a registered matrix builder by name.

    Raises KeyError if the builder name is unknown.
    """
    if name not in MATRIX_BUILDERS:
        raise KeyError(f"Unknown matrix builder: {name}")
    logger.info(
        f"Running matrix builder '{name}' with inputs={inputs} output={output_mt} params={params}"
    )
    MATRIX_BUILDERS[name](inputs, output_mt, params or {})
