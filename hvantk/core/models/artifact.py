"""Marker protocol shared by AnnotationTable / ExpressionMatrix / GeneSet.

Lets core/io dispatch without circular imports. The Artifact protocol is
intentionally tiny — every concrete artifact has a Provenance and a save()
method; everything else is type-specific.
"""
from __future__ import annotations

from pathlib import Path
from typing import Protocol, runtime_checkable

from hvantk.core.models.provenance import Provenance


@runtime_checkable
class Artifact(Protocol):
    provenance: Provenance

    def save(self, path: str | Path) -> None: ...
