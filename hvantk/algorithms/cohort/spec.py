"""Cohort manifest parsing and validation (the external-cohort contract).

A cohort manifest declares how an external cohort is presented to hvantk: which column
carries its gene identifier, which cohort-derived statistic is its prior (and in which
direction), optionally where its labels live, and optionally which of its columns form
which feature axis.

A cohort is a spine-mappable gene key PLUS a cohort-derived prior statistic (design D2).
Anything with no cohort statistic is a gene set, not a cohort. Labels are optional --
requirements are per-op, not per-cohort (D4): ``rerank`` needs labels, a credibility/veto
report does not, and extending annotation needs only the gene key.

Parsing is pure Python -- no Hail -- so manifests validate in the fast test suite. This
mirrors ``hvantk/algorithms/annotation/spec.py``, and ``CohortAxis`` deliberately mirrors
``SourceEntry``'s (axis, columns) shape so a cohort axis and a Layer-1 axis are the same
concept expressed at two layers.

Layering: this module imports stdlib + jsonschema + yaml only. It must never import
``hvantk.skills`` or ``hvantk.tools``.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import jsonschema
import yaml

LOWER_IS_BETTER = "lower_is_better"
HIGHER_IS_BETTER = "higher_is_better"

_SCHEMA_PATH = (
    Path(__file__).resolve().parents[2]
    / "resources"
    / "schemas"
    / "cohort_manifest.schema.json"
)


def _schema() -> dict:
    return json.loads(_SCHEMA_PATH.read_text())


@dataclass(frozen=True)
class CohortPrior:
    """The cohort-derived statistic that defines the ranking being re-ranked.

    ``direction`` is mandatory and never inferred: a cohort may supply a p-value
    (lower is better), a Z-score or a Bayes factor (higher is better). Guessing
    produces a silently inverted ranking -- plausible output, wrong answer.
    """

    column: str
    direction: str


@dataclass(frozen=True)
class CohortLabels:
    """Optional labels: a pointer to a GeneSetCollection JSON, never a column.

    Labels change on a different clock than burden statistics, and a collection
    records where it came from (``metadata.created_by`` / ``input_file`` /
    ``hgnc_validated``) in a way a bare column cannot.
    """

    gene_set: str
    set_name: str | None = None
    min_mapping_rate: float = 0.9


@dataclass(frozen=True)
class CohortAxis:
    """One named group of cohort columns. Mirrors ``SourceEntry``'s (axis, columns)."""

    axis: str
    columns: tuple[str, ...]


@dataclass(frozen=True)
class CohortManifest:
    """Declares how an external cohort is presented to hvantk.

    ``key`` and ``key_column`` answer two different questions -- conflating them is
    exactly the bug this pair of fields exists to prevent:

    * ``key`` is the identifier SPACE the cohort's gene key lives in --
      ``gene_id`` (Ensembl), ``hgnc_id``, or ``symbol``. This is what
      ``prepare_source`` needs to know *how* to resolve the values (which mapper
      to dispatch to); the enum stays closed to those three spine-mappable spaces.
    * ``key_column`` is the COLUMN NAME in the cohort's own table that actually
      holds those identifier values. It defaults to ``key`` when the manifest
      omits it, so a manifest whose identifier column happens to be named after
      its own id-space needs nothing extra. Real cohorts routinely don't: a
      symbol-keyed table commonly names that column ``gene``, not ``symbol`` --
      write ``key: symbol`` / ``key_column: gene`` for that case.

    ``key_column`` is always resolved to a concrete string (never ``None``) after
    construction, whether the manifest is parsed via :func:`load_cohort` or built
    directly -- no consumer has to write ``manifest.key_column or manifest.key``.
    """

    name: str
    key: str
    table: str
    prior: CohortPrior
    key_column: str | None = None
    min_mapping_rate: float = 0.9
    labels: CohortLabels | None = None
    cohort_axes: tuple[CohortAxis, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if self.key_column is None:
            object.__setattr__(self, "key_column", self.key)

    def declared_columns(self) -> tuple[str, ...]:
        """The prior column followed by every axis column, in declaration order.

        This is the exact set of columns ``attach`` selects from the cohort table;
        uniqueness across it is enforced at load time.
        """
        cols = [self.prior.column]
        for entry in self.cohort_axes:
            cols.extend(entry.columns)
        return tuple(cols)

    def axis(self, name: str) -> CohortAxis:
        for entry in self.cohort_axes:
            if entry.axis == name:
                return entry
        raise KeyError(f"no cohort axis {name!r} in manifest {self.name!r}")


def load_cohort(path: str | Path) -> CohortManifest:
    """Parse and validate a cohort manifest YAML.

    Raises
    ------
    jsonschema.ValidationError
        If the document violates the manifest schema.
    ValueError
        If two axes share a label, or if any declared column (prior or axis) is
        declared more than once.
    """
    doc = yaml.safe_load(Path(path).read_text())
    schema = _schema()
    jsonschema.validators.validator_for(schema)(schema).validate(doc)

    axes = tuple(
        CohortAxis(axis=e["axis"], columns=tuple(e["columns"]))
        for e in doc.get("cohort_axes", [])
    )
    _check_unique_axis_labels(axes)

    raw_labels = doc.get("labels")
    labels = (
        CohortLabels(
            gene_set=raw_labels["gene_set"],
            set_name=raw_labels.get("set_name"),
            min_mapping_rate=raw_labels.get("min_mapping_rate", 0.9),
        )
        if raw_labels
        else None
    )

    manifest = CohortManifest(
        name=doc["name"],
        key=doc["key"],
        table=doc["table"],
        prior=CohortPrior(
            column=doc["prior"]["column"], direction=doc["prior"]["direction"]
        ),
        key_column=doc.get("key_column", doc["key"]),
        min_mapping_rate=doc.get("min_mapping_rate", 0.9),
        labels=labels,
        cohort_axes=axes,
    )
    _check_no_duplicate_columns(manifest)
    return manifest


def _check_unique_axis_labels(axes: tuple[CohortAxis, ...]) -> None:
    seen: set[str] = set()
    for entry in axes:
        if entry.axis in seen:
            raise ValueError(
                f"duplicate cohort axis {entry.axis!r}; each axis label must be "
                "unique -- it is the handle that ablation and reporting are keyed on"
            )
        seen.add(entry.axis)


def _check_no_duplicate_columns(manifest: CohortManifest) -> None:
    """Every declared column must be declared exactly once.

    This covers both a column claimed by two axes and an axis column that collides
    with the prior column. The latter is the subtler bug: a prior statistic and a
    same-named model feature are usually different transforms of the same quantity
    (a raw p-value versus its negative log), so declaring both quietly feeds the
    untransformed value into the model.
    """
    owner: dict[str, str] = {}
    for col in (manifest.prior.column,):
        owner[col] = "prior"
    for entry in manifest.cohort_axes:
        for col in entry.columns:
            if col in owner:
                raise ValueError(
                    f"duplicate declared column {col!r}: declared by "
                    f"{owner[col]!r} and by axis {entry.axis!r}; every declared "
                    "column must appear exactly once across the prior and all "
                    "cohort axes"
                )
            owner[col] = f"axis {entry.axis}"
