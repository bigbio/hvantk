"""GeneSet: a named collection of gene identifiers with provenance."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

from hvantk.core.models.provenance import Provenance


@dataclass
class GeneSet:
    name: str
    provenance: Provenance
    _members: frozenset[str] = field(default_factory=frozenset)

    def __contains__(self, gene_id: str) -> bool:
        return gene_id in self._members

    def __len__(self) -> int:
        return len(self._members)

    def __iter__(self) -> Iterator[str]:
        return iter(self._members)

    def intersection(self, other: "GeneSet") -> "GeneSet":
        return GeneSet(
            name=f"{self.name}∩{other.name}",
            provenance=self.provenance,
            _members=self._members & other._members,
        )

    def union(self, other: "GeneSet") -> "GeneSet":
        return GeneSet(
            name=f"{self.name}∪{other.name}",
            provenance=self.provenance,
            _members=self._members | other._members,
        )

    def difference(self, other: "GeneSet") -> "GeneSet":
        return GeneSet(
            name=f"{self.name}-{other.name}",
            provenance=self.provenance,
            _members=self._members - other._members,
        )

    def to_list(self) -> list[str]:
        return list(self._members)

    def to_set(self) -> set[str]:
        return set(self._members)

    def save(self, path: str | Path) -> None:
        raise NotImplementedError("save lands in Task 15 alongside core/io")
