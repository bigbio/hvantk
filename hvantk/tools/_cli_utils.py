"""Small helpers shared by the tool CLIs (pure Python, no Hail)."""

from __future__ import annotations

import click


def parse_named_paths(values: tuple[str, ...], option: str) -> dict[str, str]:
    """Parse a repeatable ``NAME=PATH`` option into a dict.

    Raises ``click.UsageError`` naming ``option`` for a value without ``=``, an empty
    name or path, or a name given twice.
    """
    paths: dict[str, str] = {}
    for value in values:
        name, sep, path = value.partition("=")
        name, path = name.strip(), path.strip()
        if not sep or not name or not path:
            raise click.UsageError(
                f"{option} value {value!r} is malformed; expected name=path"
            )
        if name in paths:
            raise click.UsageError(
                f"{option} given twice for {name!r} "
                f"({paths[name]!r} then {path!r}); pass each name once"
            )
        paths[name] = path
    return paths
