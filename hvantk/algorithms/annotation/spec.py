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


# The `combine` vocabulary, as a dispatch table rather than a name list, so that
# SpecificitySpec's validation and reduce_matrix_to_gene's dispatch read from ONE
# definition and cannot drift into disagreement. Consumed by matrix.py; kept here
# because it is spec vocabulary, and because matrix.py already depends on this
# module's types while nothing here depends on matrix.py.
#
#   sum  -- cell-CLASS specificity (EWCE level 1): targets are subtypes of one class
#           (atrial/ventricular/Myoz2 cardiomyocytes), so pool their fractions into the
#           fraction of the gene's expression sitting in the class. A pan-class gene,
#           split across subtypes, reads high here and is missed by max.
#   mean -- average specificity across targets.
#   max  -- peak specificity to any single target.
COMBINE_REDUCERS = {
    "sum": lambda df: df.sum(axis=1),
    "mean": lambda df: df.mean(axis=1),
    "max": lambda df: df.max(axis=1),
}


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
    # Name of the emitted row-count column. The count is "source rows that survived the
    # filter", so only dbNSFP's is literally 'possible missense'; a PTM-site or eQTL-pair
    # source counts something else, and two such axes would collide on a shared name.
    # Defaults to the dbNSFP-era name so existing specs are unaffected.
    count_name: str = "n_possible_missense"


@dataclass(frozen=True)
class SpecificitySpec:
    """How a genes x groups specificity matrix becomes feature columns.

    ``emit`` controls the reduction, and defaults to the VECTOR -- one column per group.
    Reducing an atlas to a single summed scalar throws away the cross-group contrast: the
    non-target groups are computed, used as the denominator of the fraction, and discarded.
    Measured on real cohorts, that reduction cost an epilepsy axis +0.061 AUC and the
    difference between significant and not, while keeping the vector raised the
    selected-maximum null by +0.0009. So the vector is the default and a named roll-up is
    additive: give ``targets`` and you get the roll-up column IN ADDITION to the vector.

    emit
        ``"vector"`` (default) one column per group, plus the roll-up when ``targets`` is
        non-empty; ``"rollup"`` the roll-up only, which requires ``targets``.
    """

    method: str
    targets: tuple[str, ...] = ()
    combine: str = "max"
    name: str = "spec"
    emit: str = "vector"

    def __post_init__(self):
        if self.emit not in ("vector", "rollup"):
            raise ValueError(
                f"emit must be 'vector' or 'rollup'; got {self.emit!r}")
        if self.emit == "rollup" and not self.targets:
            raise ValueError(
                "emit='rollup' needs targets; with none there is nothing to roll up "
                "(an empty target set would silently sum to an all-zero column)")
        # Validated for the same reason as emit, and it matters more: the old
        # reduce_matrix_to_gene dispatched sum/mean/else-max, so an unrecognised
        # value did not raise -- it silently became 'max'. A spec author writing
        # combine: 'sm' would get peak single-target specificity where they asked
        # for the pooled class fraction: a different feature, not a degraded one.
        if self.combine not in COMBINE_REDUCERS:
            raise ValueError(
                f"combine must be one of {sorted(COMBINE_REDUCERS)}; "
                f"got {self.combine!r}"
            )


@dataclass(frozen=True)
class MatrixSpec:
    group_axis: str
    atlas: str
    tissue_tag: str | None = None
    stats: tuple[str, ...] = ()
    drop_groups: tuple[str, ...] = ()
    specificity: SpecificitySpec | None = None


@dataclass(frozen=True)
class SourceEntry:
    axis: str
    source: str
    key: str
    columns: tuple[str, ...]
    min_mapping_rate: float = 0.9
    origin: str | None = None
    ablate_separately: bool = False
    # How to reduce when several source keys map to ONE gene_id. None (default) makes
    # that an error; see prepare.COLLAPSE_REDUCERS for the vocabulary.
    collapse: str | None = None
    aggregate: AggregateSpec | None = None
    matrix: MatrixSpec | None = None


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
        count_name=raw.get("count_name", "n_possible_missense"),
    )


def _build_matrix(raw: dict | None) -> MatrixSpec | None:
    if raw is None:
        return None
    specificity_raw = raw.get("specificity")
    specificity = (
        None
        if specificity_raw is None
        else SpecificitySpec(
            method=specificity_raw["method"],
            # Optional in the schema: with no targets you get the vector alone.
            # Defaulting here rather than indexing keeps YAML able to express
            # vector-only, which is the point of the roll-up being additive.
            targets=tuple(specificity_raw.get("targets", ())),
            combine=specificity_raw.get("combine", "max"),
            name=specificity_raw.get("name", "spec"),
            # Without this, emit was Python-API-only: every YAML-driven spec
            # silently took the "vector" default and could not opt back into the
            # roll-up-only output the docstring advertises.
            emit=specificity_raw.get("emit", "vector"),
        )
    )
    return MatrixSpec(
        group_axis=raw["group_axis"],
        atlas=raw["atlas"],
        tissue_tag=raw.get("tissue_tag"),
        stats=tuple(raw.get("stats", ())),
        drop_groups=tuple(raw.get("drop_groups", ())),
        specificity=specificity,
    )


def load_spec(path: str | Path) -> FeatureSpec:
    """Read a feature-spec YAML, validate it against the schema, and return a FeatureSpec.

    Raises
    ------
    jsonschema.ValidationError
        If the document does not conform to feature_spec.schema.json.
    ValueError
        If two layer1 entries share an axis label, an entry's ``key`` and ``aggregate``
        presence disagree (``key == "variant"`` requires ``aggregate`` and vice versa), or an
        entry has a ``matrix`` block with ``key`` other than ``"symbol"`` (the reduced table
        the matrix reducer produces is symbol-keyed).
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
            collapse=e.get("collapse"),
            aggregate=_build_aggregate(e.get("aggregate")),
            matrix=_build_matrix(e.get("matrix")),
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
    for e in entries:
        if e.matrix is not None and e.key != "symbol":
            raise ValueError(
                f"entry {e.source!r} has a matrix block, which requires key: symbol "
                f"(the reduced table is symbol-keyed); got key {e.key!r}"
            )
    return FeatureSpec(name=doc["name"], layer1=entries)
