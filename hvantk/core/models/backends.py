"""
Backend declarations and algorithm metadata.

The @algorithm decorator marks functions with the backends they support
and the data format they expect, enabling the BackendRouter to select
the best execution strategy at runtime.
"""

from dataclasses import dataclass
from enum import Enum
from functools import wraps
from typing import Callable, List, Optional


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
    name : str
        Registry name for the algorithm.
    backends : list[Backend]
        Backends the algorithm supports (hard constraint for the router).
    input_format : str
        Expected input type: ``"dataframe"`` (pandas) or ``"table"`` (Hail).
    output_format : str
        Return type: ``"dataframe"`` or ``"table"``.
    key_fields : list[str] or None
        Fields to key the output Hail Table by when persisting.
    """

    backends: List[Backend]
    name: str = ""
    input_format: str = "dataframe"
    output_format: str = "dataframe"
    key_fields: Optional[List[str]] = None


_IMPLICIT_HAIL_META = AlgorithmMeta(
    backends=[Backend.HAIL],
    input_format="table",
    output_format="table",
)


def algorithm(
    backends: List[Backend],
    name: str = "",
    input_format: str = "dataframe",
    output_format: str = "dataframe",
    key_fields: Optional[List[str]] = None,
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
