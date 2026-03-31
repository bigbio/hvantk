"""
Backend declarations and algorithm metadata.

The @algorithm decorator marks functions with the backends they support
and the data format they expect, enabling the BackendRouter to select
the best execution strategy at runtime.
"""

from dataclasses import dataclass, field
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
    input_format: str = "dataframe",
    output_format: str = "dataframe",
    key_fields: Optional[List[str]] = None,
) -> Callable:
    """Declare an algorithm's supported backends and I/O formats.

    This decorator attaches an :class:`AlgorithmMeta` instance to the
    function as ``_algorithm_meta``.  It does not alter the function's
    runtime behaviour.

    Parameters
    ----------
    backends : list[Backend]
        Which backends can execute this algorithm.
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
            input_format=input_format,
            output_format=output_format,
            key_fields=key_fields,
        )

        @wraps(fn)
        def wrapper(*args, **kwargs):
            return fn(*args, **kwargs)

        wrapper._algorithm_meta = fn._algorithm_meta
        return wrapper

    return decorator


def get_algorithm_meta(fn: Callable) -> AlgorithmMeta:
    """Retrieve backend metadata from a decorated function.

    Functions without ``@algorithm`` are treated as Hail-only.
    """
    return getattr(fn, "_algorithm_meta", _IMPLICIT_HAIL_META)
