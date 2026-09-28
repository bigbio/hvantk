"""Install PEP 562 lazy exports into a package namespace.

Four hand-rolled copies of this pattern existed before #306 (``core/__init__``,
``core/models/__init__``, ``algorithms/enrichex/__init__``, ``algorithms/ptm/__init__``),
plus ``LazyGroup`` for the CLI and a bespoke ``_LazyHgcGroup``. Each one exists for the
same reason: Python imports a package before any of its submodules, so a CLI that only
needs ``package.constants`` pays for every sibling the ``__init__`` imports eagerly --
Hail (~5 s), pandas, matplotlib. This helper is the one implementation the new
conversions share; the older copies can migrate when they are next touched.
"""

from __future__ import annotations

import importlib
from typing import Callable, Mapping


def install_lazy_exports(
    namespace: dict,
    exports: Mapping[str, str],
    *,
    missing_hint: Callable[[ModuleNotFoundError], "str | None"] | None = None,
) -> None:
    """Make ``exports`` (attribute name -> defining module) resolve on first access.

    Adds ``__getattr__`` and ``__dir__`` to ``namespace`` (a package's ``globals()``).
    A resolved value is cached in the namespace, so the import cost is paid once. When
    the defining module cannot be imported because a dependency is absent,
    ``missing_hint(exc)`` may return the message of the ``ImportError`` to raise instead
    (naming the extra to install); returning ``None`` propagates the original error.
    Attribute names that are not exports but name a submodule of the package are
    imported and returned, so ``import pkg as p; p.submodule`` keeps working for a
    submodule nobody imported yet.
    """
    package = namespace["__name__"]
    exports = dict(exports)

    def __getattr__(name: str):
        module_name = exports.get(name)
        if module_name is None:
            qualified = f"{package}.{name}"
            try:
                return importlib.import_module(qualified)
            except ModuleNotFoundError as exc:
                # exc.name is the qualified submodule itself when the parent package
                # imports fine but has no such submodule; it is a *prefix* of qualified
                # when the parent package cannot be found at all (import machinery fails
                # on the parent before it ever gets to the submodule). Both mean "there is
                # no such submodule" from this package's point of view; anything else is
                # an unrelated import failure nested inside a submodule that does exist,
                # and must propagate rather than be mistaken for a missing attribute.
                if exc.name == qualified or qualified.startswith(f"{exc.name}."):
                    raise AttributeError(f"module {package!r} has no attribute {name!r}") from None
                raise
        try:
            value = getattr(importlib.import_module(module_name), name)
        except ModuleNotFoundError as exc:
            hint = missing_hint(exc) if missing_hint is not None else None
            if hint:
                raise ImportError(hint) from exc
            raise
        namespace[name] = value
        return value

    def __dir__():
        return sorted(set(namespace) | set(exports))

    namespace["__getattr__"] = __getattr__
    namespace["__dir__"] = __dir__
