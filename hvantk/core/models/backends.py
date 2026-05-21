"""
Backend declarations and algorithm metadata.

The @algorithm decorator marks functions with the backends they support
and the data format they expect, enabling the BackendRouter to select
the best execution strategy at runtime.
"""

from dataclasses import dataclass, field
from enum import Enum
from functools import wraps
from typing import Callable


class Backend(Enum):
    """Supported compute/IO backends."""

    HAIL = "hail"
    PANDAS = "pandas"
    DUCKDB = "duckdb"


@dataclass
class AlgorithmMeta:
    """Metadata attached to algorithm functions by the @algorithm decorator.

    Attributes
    ----------
    backends : list[Backend]
        Backends the algorithm supports (hard constraint for the router).
    name : str
        Registry name. Defaults to the decorated function's name.
    inputs : dict[str, type]
        Typed input schema (param_name -> artifact type). Empty by default.
        Populated by the decorator's ``inputs=`` argument.
    outputs : dict[str, type]
        Typed output schema (name -> artifact type). Empty by default.
    required_backend : str | None
        If set, declares a hard backend requirement (escape hatch for
        algorithms that can't be backend-agnostic, e.g. genotype workflows
        consuming raw hl.MatrixTable). Use sparingly; document the
        justification in the algorithm's docstring.
    input_format : str
        Legacy: ``"dataframe"`` or ``"table"``. Retained for backward compat.
    output_format : str
        Legacy: ``"dataframe"`` or ``"table"``. Retained for backward compat.
    key_fields : list[str] | None
        Legacy: Hail Table key fields for persisting the output.
    """

    backends: list[Backend]
    name: str = ""
    inputs: dict[str, type] = field(default_factory=dict)
    outputs: dict[str, type] = field(default_factory=dict)
    required_backend: str | None = None
    input_format: str = "dataframe"
    output_format: str = "dataframe"
    key_fields: list[str] | None = None


_IMPLICIT_HAIL_META = AlgorithmMeta(
    backends=[Backend.HAIL],
    inputs={},
    outputs={},
    required_backend=None,
    input_format="table",
    output_format="table",
)


def algorithm(
    backends: list[Backend],
    name: str = "",
    inputs: dict[str, type] | None = None,
    outputs: dict[str, type] | None = None,
    required_backend: str | None = None,
    input_format: str = "dataframe",
    output_format: str = "dataframe",
    key_fields: list[str] | None = None,
) -> Callable:
    """Declare an algorithm's supported backends and I/O formats.

    This decorator attaches an :class:`AlgorithmMeta` instance to the
    function as ``_algorithm_meta``.

    Provenance chaining
    -------------------
    When the decorated function is called with Artifact inputs (objects
    with a ``.provenance`` attribute that is a ``Provenance`` instance) and
    returns an Artifact output, the output's ``provenance.parents`` is
    automatically set to a tuple of the input provenances. This preserves
    the build graph for downstream drift detection and reproducibility.
    If the algorithm already stamped non-empty parents on the output, the
    existing chain is preserved (the decorator does not overwrite it).

    Parameters
    ----------
    backends : list[Backend]
        Which backends can execute this algorithm.
    name : str, optional
        Registry name for the algorithm. Defaults to the function name.
    inputs : dict[str, type], optional
        Typed input schema (param_name -> artifact type).
    outputs : dict[str, type], optional
        Typed output schema (name -> artifact type).
    required_backend : str | None, optional
        Hard backend requirement; use sparingly.
    input_format : str
        ``"dataframe"`` or ``"table"`` — what the function receives.
    output_format : str
        ``"dataframe"`` or ``"table"`` — what the function returns.
    key_fields : list[str], optional
        Hail Table key fields for persisting the output.
    """

    def decorator(fn: Callable) -> Callable:
        fn._algorithm_meta = AlgorithmMeta(
            backends=backends,
            name=name or fn.__name__,
            inputs=inputs or {},
            outputs=outputs or {},
            required_backend=required_backend,
            input_format=input_format,
            output_format=output_format,
            key_fields=key_fields,
        )

        @wraps(fn)
        def wrapper(*args, **kwargs):
            # Collect input artifact provenances by introspecting args/kwargs.
            # We treat anything with a .provenance attribute that's a Provenance
            # instance as an Artifact. This covers AnnotationTable, ExpressionMatrix,
            # GeneSet, and any future artifact types without requiring an explicit
            # isinstance check (which would force an import cycle).
            from hvantk.core.models.provenance import Provenance

            input_provs: list = []
            for arg in args:
                prov = getattr(arg, "provenance", None)
                if isinstance(prov, Provenance):
                    input_provs.append(prov)
            for kwarg in kwargs.values():
                prov = getattr(kwarg, "provenance", None)
                if isinstance(prov, Provenance):
                    input_provs.append(prov)

            result = fn(*args, **kwargs)

            # If the result is an artifact-like object and didn't already chain
            # provenance internally, attach the inputs' provenances as its parents.
            if input_provs and hasattr(result, "provenance"):
                result_prov = getattr(result, "provenance", None)
                if isinstance(result_prov, Provenance) and result_prov.parents == ():
                    import dataclasses
                    result.provenance = dataclasses.replace(
                        result_prov, parents=tuple(input_provs)
                    )

            return result

        wrapper._algorithm_meta = fn._algorithm_meta
        return wrapper

    return decorator


def get_algorithm_meta(fn: Callable) -> AlgorithmMeta:
    """Retrieve backend metadata from a decorated function.

    Functions without ``@algorithm`` are treated as Hail-only.
    """
    return getattr(fn, "_algorithm_meta", _IMPLICIT_HAIL_META)
