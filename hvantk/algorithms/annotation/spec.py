"""Feature-spec parsing and validation.

A feature spec declares, per source, which columns become gene-level features and how the
source's key maps onto the spine's ``gene_id``. For an ML training set (``hvantk
training-set build``) the same file also declares the ``label``, ``variant``, ``gene`` and
``gene_lists`` sections and a ``tier`` per layer1 axis. A ``gene_id``-keyed entry joins directly;
other keys (e.g. ``symbol``) are mapped onto ``gene_id``, and ``variant``-keyed entries are
aggregated through an ``aggregate`` block. Parsing is pure Python -- no Hail -- so specs
validate in the fast test suite.
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
    That contrast can carry the signal: in measurement, the reduction lost real AUC, while
    keeping the vector barely moved the selected-maximum null. So the vector is the default
    and a named roll-up is additive: give ``targets`` and you get the roll-up column IN
    ADDITION to the vector.

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
            raise ValueError(f"emit must be 'vector' or 'rollup'; got {self.emit!r}")
        if self.emit == "rollup" and not self.targets:
            raise ValueError(
                "emit='rollup' needs targets; with none there is nothing to roll up "
                "(an empty target set would silently sum to an all-zero column)"
            )
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
    # Leakage tier of the axis's columns in an ML training set (see
    # hvantk.algorithms.annotation.tiers); unused by annotate prepare/compose.
    tier: str | None = None


# Columns every training set carries ahead of its features (never features themselves).
METADATA_COLUMNS = (
    "vid",
    "locus",
    "alleles",
    "gene_id",
    "gene_symbol",
    "clinvar_clnsig",
    "clinvar_review_status",
    "clinvar_review_stars",
    "consequence",
    "label",
)


@dataclass(frozen=True)
class LabelSpec:
    """Where a training set's rows and 1/0 label come from (``label:``)."""

    input: str
    source: str
    provenance: tuple[str, ...]
    pathogenic: tuple[str, ...]
    benign: tuple[str, ...]
    min_review_stars: int = 1
    consequence: str = "missense_variant"
    snv_only: bool = True


@dataclass(frozen=True)
class VariantColumn:
    """One variant-level output column.

    For an ``allele`` join, ``field`` (or ``fields`` reduced by ``reduce``) is a dotted
    path into the source table; ``cast: float64`` parses a string field. An ``interval``
    or ``codon`` join yields a boolean flag and names no field.
    """

    name: str
    field: str | None = None
    fields: tuple[str, ...] = ()
    cast: str | None = None
    reduce: str | None = None
    tier: str | None = None
    description: str = ""


@dataclass(frozen=True)
class VariantSource:
    """One variant-level source joined onto the labelled rows (an item of ``variant:``).

    ``join`` is ``allele`` (same ``locus``/``alleles``), ``interval`` (the locus lies in
    one of the source's intervals) or ``codon`` (the locus lies in a site's codon).
    ``gene_field`` names the source's ``;``-separated Ensembl gene ids, used to choose
    the gene of a variant ClinVar assigns to several genes.
    """

    input: str
    source: str
    join: str
    columns: tuple[VariantColumn, ...]
    gene_field: str | None = None


@dataclass(frozen=True)
class SpineColumn:
    name: str
    tier: str
    description: str = ""


@dataclass(frozen=True)
class GeneSpec:
    """How variants reach genes (``gene:``): the HGNC input and the spine columns kept."""

    hgnc_input: str
    spine_columns: tuple[SpineColumn, ...] = ()


@dataclass(frozen=True)
class GeneList:
    column: str
    description: str = ""


@dataclass(frozen=True)
class FeatureSpec:
    name: str
    layer1: tuple[SourceEntry, ...] = field(default_factory=tuple)
    label: LabelSpec | None = None
    variant: tuple[VariantSource, ...] = ()
    gene: GeneSpec | None = None
    gene_lists: tuple[GeneList, ...] = ()

    def entry(self, axis: str) -> SourceEntry:
        for e in self.layer1:
            if e.axis == axis:
                return e
        raise KeyError(f"no layer1 entry for axis {axis!r}")

    def input_names(self) -> tuple[str, ...]:
        """The ``training-set build --input`` names this spec needs, in spec order."""
        names = [self.label.input] if self.label is not None else []
        names += [v.input for v in self.variant]
        if self.gene is not None:
            names.append(self.gene.hgnc_input)
        return tuple(names)


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


def _build_label(raw: dict | None) -> LabelSpec | None:
    if raw is None:
        return None
    return LabelSpec(
        input=raw["input"],
        source=raw["source"],
        provenance=tuple(raw["provenance"]),
        pathogenic=tuple(raw["pathogenic"]),
        benign=tuple(raw["benign"]),
        min_review_stars=raw.get("min_review_stars", 1),
        consequence=raw.get("consequence", "missense_variant"),
        snv_only=raw.get("snv_only", True),
    )


def _build_variant_source(raw: dict) -> VariantSource:
    columns = tuple(
        VariantColumn(
            name=name,
            field=col.get("field"),
            fields=tuple(col.get("fields", ())),
            cast=col.get("cast"),
            reduce=col.get("reduce"),
            tier=col.get("tier"),
            description=col.get("description", ""),
        )
        for name, col in raw["columns"].items()
    )
    source = VariantSource(
        input=raw["input"],
        source=raw["source"],
        join=raw["join"],
        columns=columns,
        gene_field=raw.get("gene_field"),
    )
    for col in columns:
        where = f"variant source {source.input!r} column {col.name!r}"
        if source.join == "allele":
            if (col.field is None) == (not col.fields):
                raise ValueError(
                    f"{where}: an allele join needs exactly one of field/fields"
                )
        elif col.field is not None or col.fields or col.cast or col.reduce:
            raise ValueError(
                f"{where}: the {source.join} join yields a flag; field, fields, cast "
                "and reduce do not apply"
            )
        if bool(col.fields) != (col.reduce is not None):
            raise ValueError(f"{where}: fields needs reduce, and reduce needs fields")
    if source.join != "allele" and len(columns) != 1:
        raise ValueError(
            f"variant source {source.input!r}: the {source.join} join yields one flag; "
            f"declare exactly one column (got {len(columns)})"
        )
    if source.gene_field is not None and source.join != "allele":
        raise ValueError(
            f"variant source {source.input!r}: gene_field needs an allele join"
        )
    return source


def _build_gene(raw: dict | None) -> GeneSpec | None:
    if raw is None:
        return None
    return GeneSpec(
        hgnc_input=raw["hgnc_input"],
        spine_columns=tuple(
            SpineColumn(name=name, tier=c["tier"], description=c.get("description", ""))
            for name, c in (raw.get("spine_columns") or {}).items()
        ),
    )


def _check_training_sections(spec: FeatureSpec) -> None:
    """Cross-section checks of the training-set sections (pure Python)."""
    inputs = spec.input_names()
    repeated = sorted({n for n in inputs if inputs.count(n) > 1})
    if repeated:
        raise ValueError(
            "input name(s) used twice in the spec (each is one --input): "
            + ", ".join(repeated)
        )
    with_gene = [v.input for v in spec.variant if v.gene_field is not None]
    if len(with_gene) > 1:
        raise ValueError(
            "gene_field is declared on several variant sources "
            f"({', '.join(with_gene)}); declare it on one"
        )
    owner: dict[str, str] = {c: "metadata" for c in METADATA_COLUMNS}
    named = [(g.column, "gene_lists") for g in spec.gene_lists]
    named += [
        (c.name, f"variant source {v.input!r}") for v in spec.variant for c in v.columns
    ]
    if spec.gene is not None:
        named += [(c.name, "gene.spine_columns") for c in spec.gene.spine_columns]
    for column, where in named:
        if column in owner:
            raise ValueError(
                f"column {column!r} is declared by {owner[column]} and by {where}; "
                "training-set column names must be unique"
            )
        owner[column] = where


def load_spec(path: str | Path) -> FeatureSpec:
    """Read a feature-spec YAML, validate it against the schema, and return a FeatureSpec.

    Raises
    ------
    jsonschema.ValidationError
        If the document does not conform to feature_spec.schema.json.
    ValueError
        If two layer1 entries share an axis label, an entry's ``key`` and ``aggregate``
        presence disagree (``key == "variant"`` requires ``aggregate`` and vice versa), or an
        entry has a ``matrix`` block with ``key`` other than ``"symbol"`` or ``"gene_id"``
        (the reduced table is keyed by the atlas's gene identifiers).
        Also, for the training-set sections: an allele-join column without exactly one of
        ``field``/``fields``, a flag join naming a field, ``fields`` without ``reduce`` (or
        the reverse), ``gene_field`` on a non-allele join or on several sources, an input
        name used twice, or a column name used twice (metadata columns included).
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
            tier=e.get("tier"),
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
        if e.matrix is not None and e.key not in ("symbol", "gene_id"):
            raise ValueError(
                f"entry {e.source!r} has a matrix block, which requires key: symbol "
                f"or key: gene_id (the atlas's genes are keyed by symbol or by "
                f"Ensembl gene id); got key {e.key!r}"
            )
    # Matrix vector columns are named {atlas}_{sanitized_group} with no axis component
    # (matrix.reduce_matrix_to_gene), so two matrix entries over the SAME atlas emit the
    # same column names for every group they share. Since emit="vector" is the default,
    # that is now the common case -- e.g. splitting one single-cell atlas into a
    # cardiomyocyte axis and a fibroblast axis. compose() would catch it, but only after
    # both sources had been read, reduced and mapped onto the spine; here it costs
    # nothing and names the fix.
    atlases: dict = {}
    for e in entries:
        if e.matrix is None:
            continue
        if e.matrix.atlas in atlases:
            raise ValueError(
                f"axes {atlases[e.matrix.atlas]!r} and {e.axis!r} both declare matrix "
                f"atlas {e.matrix.atlas!r}; their per-group specificity columns would "
                "collide, since a vector column is named {atlas}_{group}. Give each "
                "matrix axis a distinct atlas label, or set emit: rollup on all but one"
            )
        atlases[e.matrix.atlas] = e.axis
    spec = FeatureSpec(
        name=doc["name"],
        layer1=entries,
        label=_build_label(doc.get("label")),
        variant=tuple(_build_variant_source(v) for v in doc.get("variant", ())),
        gene=_build_gene(doc.get("gene")),
        gene_lists=tuple(
            GeneList(column=name, description=(g or {}).get("description", ""))
            for name, g in (doc.get("gene_lists") or {}).items()
        ),
    )
    _check_training_sections(spec)
    return spec
