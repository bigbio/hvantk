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
class SourceEntry:
    axis: str
    source: str
    key: str
    columns: tuple[str, ...]
    min_mapping_rate: float = 0.9
    origin: str | None = None
    ablate_separately: bool = False


@dataclass(frozen=True)
class FeatureSpec:
    name: str
    layer1: tuple[SourceEntry, ...] = field(default_factory=tuple)

    def entry(self, axis: str) -> SourceEntry:
        for e in self.layer1:
            if e.axis == axis:
                return e
        raise KeyError(f"no layer1 entry for axis {axis!r}")


def load_spec(path: str | Path) -> FeatureSpec:
    """Read a feature-spec YAML, validate it against the schema, and return a FeatureSpec.

    Raises
    ------
    jsonschema.ValidationError
        If the document does not conform to feature_spec.schema.json.
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
        )
        for e in doc["layer1"]
    )
    return FeatureSpec(name=doc["name"], layer1=entries)
