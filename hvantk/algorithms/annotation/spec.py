"""Feature-spec parsing and validation.

A feature spec declares, per source, which columns become gene-level features and how the
source's key maps onto the spine's ``gene_id``. P2a supports only ``gene_id``-keyed
entries (direct); P2c widens the key types and adds aggregation transforms. Parsing is
pure Python -- no Hail -- so specs validate in the fast test suite.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import jsonschema
import yaml

_SCHEMA_PATH = (
    Path(__file__).resolve().parents[2]
    / "resources"
    / "schemas"
    / "feature_spec.schema.json"
)


def _schema() -> dict:
    return json.loads(_SCHEMA_PATH.read_text())


@dataclass(frozen=True)
class ScoreSpec:
    name: str
    column: str
    stats: tuple[str, ...]


@dataclass(frozen=True)
class AggregateSpec:
    by: str
    to: str
    scores: tuple[ScoreSpec, ...]
    filter: str | None = None
    reduce: str = "max"


@dataclass(frozen=True)
class SourceEntry:
    axis: str
    source: str
    key: str
    columns: tuple[str, ...]
    min_mapping_rate: float = 0.9
    origin: str | None = None
    ablate_separately: bool = False
    aggregate: AggregateSpec | None = None


@dataclass(frozen=True)
class FeatureSpec:
    name: str
    layer1: tuple[SourceEntry, ...] = field(default_factory=tuple)

    def entry(self, axis: str) -> SourceEntry:
        for e in self.layer1:
            if e.axis == axis:
                return e
        raise KeyError(f"no layer1 entry for axis {axis!r}")


def _build_aggregate(raw: dict | None) -> AggregateSpec | None:
    if raw is None:
        return None
    scores = tuple(
        ScoreSpec(name=name, column=sc["column"], stats=tuple(sc["stats"]))
        for name, sc in raw["scores"].items()
    )
    return AggregateSpec(
        by=raw["by"],
        to=raw["to"],
        scores=scores,
        filter=raw.get("filter"),
        reduce=raw.get("reduce", "max"),
    )


def load_spec(path: str | Path) -> FeatureSpec:
    """Read a feature-spec YAML, validate it against the schema, and return a FeatureSpec.

    Raises
    ------
    jsonschema.ValidationError
        If the document does not conform to feature_spec.schema.json.
    ValueError
        If two layer1 entries share an axis label, or an entry's ``key`` and ``aggregate``
        presence disagree (``key == "variant"`` requires ``aggregate`` and vice versa).
    """
    doc = yaml.safe_load(Path(path).read_text())
    jsonschema.validate(doc, _schema())  # raises ValidationError on failure
    entries = tuple(
        SourceEntry(
            axis=e["axis"],
            source=e["source"],
            key=e["key"],
            columns=tuple(e["columns"]),
            min_mapping_rate=e.get("min_mapping_rate", 0.9),
            origin=e.get("origin"),
            ablate_separately=e.get("ablate_separately", False),
            aggregate=_build_aggregate(e.get("aggregate")),
        )
        for e in doc["layer1"]
    )
    axes = [e.axis for e in entries]
    duplicates = sorted({a for a in axes if axes.count(a) > 1})
    if duplicates:
        raise ValueError(
            "duplicate axis label(s) in spec (each layer1 axis must be unique, it is the "
            f"CLI addressing key): {', '.join(duplicates)}"
        )
    for e in entries:
        if (e.key == "variant") != (e.aggregate is not None):
            raise ValueError(
                f"entry {e.source!r}: key 'variant' requires an aggregate block and vice versa"
            )
    return FeatureSpec(name=doc["name"], layer1=entries)
