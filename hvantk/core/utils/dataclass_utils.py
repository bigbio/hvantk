"""Rebuild dataclass configs from saved dictionaries."""

from dataclasses import fields
from typing import Any, Mapping, Type, TypeVar

T = TypeVar("T")


def dataclass_from_dict(cls: Type[T], data: Mapping[str, Any]) -> T:
    """Build ``cls`` from ``data``, ignoring keys ``cls`` does not define.

    Pipeline state files store ``asdict(config)``. Dropping the keys a newer
    version no longer has lets a state file written by an older version still load.
    """
    known = {f.name for f in fields(cls)}
    return cls(**{k: v for k, v in data.items() if k in known})
