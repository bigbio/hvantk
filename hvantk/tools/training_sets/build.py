"""Assemble an ML training dataset from a feature spec (``hvantk training-set build``).

Rows are ClinVar variants labelled 1/0 (``label:``); variant-level columns come from the
``variant:`` sources (allele-keyed scores and frequencies, interval and codon flags); each
row reaches a gene through ClinVar's GENEINFO (``gene:``) and takes that gene's columns
from the ``annotate compose`` matrix; ``gene_lists:`` add membership flags. Every feature
column gets a leakage tier (:mod:`hvantk.algorithms.annotation.tiers`); a column whose
tier cannot be derived is left out and listed in the manifest.

Outputs (:func:`write_training_set`): ``<prefix>.ht``, a deterministic ``<prefix>.tsv.gz``,
``<prefix>.dictionary.tsv`` (one row per column) and ``<prefix>.manifest.json`` (counts,
tiers, exclusions, gene-list mapping, input provenance).

Layering: tools/ composes the clinvar skill (labels) with algorithms/ (spec, tiers,
variant transforms, gene mapping).
"""

from __future__ import annotations

import dataclasses
import importlib.metadata
import json
import logging
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import hail as hl
import pandas as pd

from hvantk.algorithms.annotation.mapping import GeneIdMapper
from hvantk.algorithms.annotation.spec import METADATA_COLUMNS, FeatureSpec
from hvantk.algorithms.annotation.tiers import TierInput, resolve_tiers
from hvantk.algorithms.annotation.transforms import output_name
from hvantk.algorithms.annotation.variant import (
    choose_candidate,
    codon_sites,
    column_expr,
    field_expr,
)
from hvantk.core.models.provenance import Provenance
from hvantk.skills.clinvar.streamers import label_clinvar_variants

logger = logging.getLogger(__name__)

SCHEMA_ID = "training-set-v1"
# Label sources this builder knows how to label (label_clinvar_variants).
LABEL_SOURCES = ("clinvar:variants",)
DICTIONARY_COLUMNS = (
    "column",
    "legacy_name",
    "training_feature",
    "tier",
    "level",
    "source",
    "source_version",
    "description",
    "non_null_rate",
    "true_rate",
)
# Metadata columns after the (locus, alleles) key, in output order.
_ROW_METADATA = tuple(c for c in METADATA_COLUMNS if c not in ("locus", "alleles"))
_PRESENT_TIER = "L1"


@dataclass(frozen=True)
class TrainingSetResult:
    """A built training set, not yet written: see :func:`write_training_set`."""

    table: hl.Table
    dictionary: pd.DataFrame
    manifest: dict
    provenance: Provenance
    work_path: str


@dataclass(frozen=True)
class _Column:
    """Everything the dictionary and the tier rule need to know about one column."""

    name: str
    origin: str  # metadata | gene_list | variant | matrix (spine and layer1 columns)
    level: str
    source: str
    description: str
    feature: bool = True
    input: str | None = None
    trained_on: frozenset[str] | None = None
    declared: str | None = None


def build_training_set(
    spec: FeatureSpec,
    inputs: Mapping[str, hl.Table],
    gene_matrix: hl.Table,
    gene_axes: Mapping[str, Sequence[str]],
    gene_lists: Mapping[str, Sequence[str]],
    hgnc: Any,
    score_provenance: Mapping[str, Mapping[str, Sequence[str]]],
    input_provenance: Mapping[str, Provenance | None],
    *,
    checkpoint: str,
    spec_sha256: str = "",
    legacy_names: Mapping[str, str] | None = None,
) -> TrainingSetResult:
    """Build the training set ``spec`` declares.

    Parameters
    ----------
    spec
        A feature spec with ``label`` and ``gene`` sections.
    inputs
        One Hail table per ``spec.input_names()``.
    gene_matrix, gene_axes
        The ``annotate compose`` output (keyed ``gene_id``) and, from its manifest, the
        columns each ``layer1`` axis contributed.
    gene_lists
        Gene symbols per ``gene_lists`` column.
    hgnc
        An ``HGNCGeneCatalogStreamer`` over the ``gene.hgnc_input`` table.
    score_provenance
        Per variant ``source``, its plugin's ``scores`` declarations (field -> trained_on).
    input_provenance
        Per input name, the input's provenance (``None`` when it has none).
    checkpoint
        Where the intermediate table is written (removed by :func:`write_training_set`).
    spec_sha256
        Hash of the spec file, recorded in the manifest and the output's provenance.
    legacy_names
        Optional column -> legacy name map for the dictionary's ``legacy_name``.

    Raises
    ------
    ValueError
        When the spec, inputs, gene lists or gene matrix disagree; checked before any
        Hail job runs.
    KeyError
        When a spec field path does not exist in its source table.
    """
    check_inputs(spec, inputs, gene_matrix, gene_axes, gene_lists)
    columns = _plan_columns(spec, gene_axes, score_provenance, input_provenance)
    tiers = resolve_tiers(
        [TierInput(c.name, c.trained_on, c.declared) for c in columns if c.feature],
        spec.label.provenance,
    )
    excluded = [
        {
            "column": c.name,
            "source": c.source,
            "reason": "no declared trained_on and no tier in the spec",
        }
        for c in columns
        if c.feature and tiers[c.name] is None
    ]
    kept = [c for c in columns if not c.feature or tiers[c.name] is not None]
    kept_names = {c.name for c in kept}

    label = spec.label
    raw = inputs[label.input]
    labelled = label_clinvar_variants(
        raw,
        pathogenic=label.pathogenic,
        benign=label.benign,
        consequence=label.consequence,
        snv_only=label.snv_only,
    )
    is_kept = hl.is_defined(labelled.label) & (
        labelled.clinvar_review_stars >= label.min_review_stars
    )
    funnel = labelled.aggregate(
        hl.struct(
            after_variant_filter=hl.agg.count(),
            labelled=hl.agg.count_where(hl.is_defined(labelled.label)),
            labelled_min_stars=hl.agg.count_where(is_kept),
        )
    )
    counts: dict = {"clinvar_rows": raw.count(), **dict(funnel)}
    ht = labelled.filter(is_kept)

    # --- variant sources ---------------------------------------------------------------
    preferring = _preferring_source(spec)
    allele_inputs = []
    for v in spec.variant:
        src = inputs[v.input]
        names = [c.name for c in v.columns if c.name in kept_names]
        if v.join == "allele":
            exprs = {c.name: column_expr(src, c) for c in v.columns if c.name in names}
            if v.gene_field is not None:
                exprs["_preferred_genes"] = hl.set(
                    hl.or_else(field_expr(src, v.gene_field), "")
                    .split(";")
                    .filter(lambda g: (g != "") & (g != "."))
                )
            joined = src.select(**exprs)[ht.key]
            ht = ht.annotate(
                **{name: joined[name] for name in exprs},
                **{f"_in_{v.input}": hl.is_defined(joined)},
            )
            allele_inputs.append(v.input)
        elif v.join == "interval":
            ht = ht.annotate(**{n: hl.is_defined(src[ht.locus]) for n in names})
        else:  # codon
            sites = codon_sites(src)
            ht = ht.annotate(**{n: hl.is_defined(sites[ht.locus]) for n in names})

    # --- gene choice ---------------------------------------------------------------------
    hgnc_ht = inputs[spec.gene.hgnc_input]
    entrez_rows = (
        hgnc_ht.filter(
            hl.is_defined(hgnc_ht.entrez_id)
            & (hgnc_ht.entrez_id != "")
            & hl.is_defined(hgnc_ht.ensembl_gene_id)
            & (hgnc_ht.ensembl_gene_id != "")
        )
        .key_by()
        .select("entrez_id", "ensembl_gene_id")
        .collect()
    )
    entrez_to_gene: dict[str, str] = {}
    for row in sorted(entrez_rows, key=lambda r: (r.entrez_id, r.ensembl_gene_id)):
        entrez_to_gene.setdefault(row.entrez_id, row.ensembl_gene_id)
    lookup = hl.literal(entrez_to_gene, hl.tdict(hl.tstr, hl.tstr))
    ht = ht.annotate(_candidates=ht.geneinfo_entrez.map(lambda e: lookup.get(e)))
    preferred = (
        ht._preferred_genes if preferring is not None else hl.missing(hl.tset(hl.tstr))
    )
    ht = ht.annotate(_choice=choose_candidate(ht._candidates, preferred))
    index = ht._choice.index
    ht = ht.annotate(
        gene_id=hl.or_missing(hl.is_defined(index), ht._candidates[index]),
        gene_symbol=hl.or_missing(hl.is_defined(index), ht.geneinfo_symbols[index]),
        _gene_rule=ht._choice.rule,
    )

    # --- gene matrix and gene lists ------------------------------------------------------
    matrix_columns = [c.name for c in kept if c.origin == "matrix"]
    joined = gene_matrix.select(*matrix_columns)[ht.gene_id]
    ht = ht.annotate(**{name: joined[name] for name in matrix_columns})

    mapper = GeneIdMapper(hgnc, gene_matrix.gene_id.collect())
    gene_list_report = {}
    for gl in spec.gene_lists:
        symbols = list(gene_lists[gl.column])
        mapping, _ = mapper.from_symbols(symbols, source=gl.column)
        gene_list_report[gl.column] = _gene_list_report(symbols, mapping)
        ids = {g for g in mapping.values() if g}
        members = hl.literal(ids, hl.tset(hl.tstr))
        ht = ht.annotate(**{gl.column: hl.or_else(members.contains(ht.gene_id), False)})

    # --- select, checkpoint, measure -----------------------------------------------------
    out = [c.name for c in kept]
    internal = [f"_in_{i}" for i in allele_inputs] + ["_gene_rule"]
    work = ht.select(*out, *internal).checkpoint(checkpoint, overwrite=True)

    all_columns = ["locus", "alleles", *out]
    booleans = [c for c in out if work[c].dtype == hl.tbool]
    gene_list_columns = [gl.column for gl in spec.gene_lists]
    agg = work.aggregate(
        hl.struct(
            n=hl.agg.count(),
            label=hl.agg.counter(work.label),
            label_by_gene_list=hl.struct(
                **{
                    g: hl.agg.filter(work[g], hl.agg.counter(work.label))
                    for g in gene_list_columns
                }
            ),
            review_stars=hl.agg.counter(work.clinvar_review_stars),
            in_input=hl.struct(
                **{i: hl.agg.count_where(work[f"_in_{i}"]) for i in allele_inputs}
            ),
            no_gene_id=hl.agg.count_where(hl.is_missing(work.gene_id)),
            gene_rule=hl.agg.counter(work._gene_rule),
            non_null=hl.struct(
                **{c: hl.agg.count_where(hl.is_defined(work[c])) for c in all_columns}
            ),
            true=hl.struct(
                **{c: hl.agg.count_where(hl.or_else(work[c], False)) for c in booleans}
            ),
        )
    )
    n = agg.n
    counts.update(
        rows=n,
        label=_str_keys(agg.label),
        label_by_gene_list={
            g: _str_keys(agg.label_by_gene_list[g]) for g in gene_list_columns
        },
        review_stars=_str_keys(agg.review_stars),
        **{f"in_{i}": agg.in_input[i] for i in allele_inputs},
        no_gene_id=agg.no_gene_id,
        gene_rule=_str_keys(agg.gene_rule),
    )

    legacy_names = legacy_names or {}
    by_name = {c.name: c for c in kept}
    by_name.update({k: _metadata_column(k, spec) for k in ("locus", "alleles")})
    rows = []
    for name in all_columns:
        c = by_name[name]
        rows.append(
            {
                "column": name,
                "legacy_name": legacy_names.get(name, ""),
                "training_feature": "yes" if c.feature else "no",
                "tier": (tiers[name] or "") if c.feature else "",
                "level": c.level,
                "source": c.source,
                "source_version": _source_version(input_provenance, c.input),
                "description": c.description,
                "non_null_rate": round(agg.non_null[name] / n, 4) if n else 0.0,
                "true_rate": (round(agg.true[name] / n, 4) if n else 0.0)
                if name in booleans
                else "",
            }
        )
    dictionary = pd.DataFrame(rows, columns=list(DICTIONARY_COLUMNS))

    manifest = {
        "spec": {"name": spec.name, "sha256": spec_sha256},
        "counts": counts,
        "tiers": {
            t: sum(1 for c in kept if c.feature and tiers[c.name] == t)
            for t in ("L1", "L2", "L3", "L4")
        },
        "excluded": excluded,
        "gene_choice": {
            "rule": "first GENEINFO gene that the preferring source also assigns, "
            "else the first GENEINFO gene",
            "preferring_source": None
            if preferring is None
            else f"{preferring.source} ({preferring.gene_field})",
        },
        "gene_lists": gene_list_report,
        "input_provenance": {
            name: None if p is None else _jsonable(dataclasses.asdict(p))
            for name, p in input_provenance.items()
        },
    }
    # The code that built it, found the way every plugin build finds it.
    from hvantk.core.plugin.run_builder import resolve_builder_commit

    provenance = Provenance(
        plugin="hvantk",
        dataset=f"training-set:{spec.name}",
        plugin_version=_hvantk_version(),
        source_fingerprint=f"spec-sha256:{spec_sha256}",
        schema_id=SCHEMA_ID,
        build_timestamp=datetime.now(timezone.utc),
        builder_commit=resolve_builder_commit(),
        parents=tuple(p for p in input_provenance.values() if p is not None),
        build_parameters={"spec": spec.name},
    )
    return TrainingSetResult(
        table=work.select(*out),
        dictionary=dictionary,
        manifest=manifest,
        provenance=provenance,
        work_path=checkpoint,
    )


def write_training_set(result: TrainingSetResult, prefix: str) -> dict[str, str]:
    """Write ``<prefix>.ht``, ``.tsv.gz``, ``.dictionary.tsv`` and ``.manifest.json``.

    The ``.ht`` and ``.tsv.gz`` carry provenance sidecars; the TSV is the deterministic
    export of :func:`hvantk.core.io.save`. The intermediate table is removed afterwards.
    Returns the written paths by kind.
    """
    from hvantk.core.io import save
    from hvantk.core.models.annotation_table import AnnotationTable

    paths = {
        "table": f"{prefix}.ht",
        "tsv": f"{prefix}.tsv.gz",
        "dictionary": f"{prefix}.dictionary.tsv",
        "manifest": f"{prefix}.manifest.json",
    }
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    save(
        AnnotationTable.from_hail(result.table, provenance=result.provenance),
        paths["table"],
    )
    save(
        AnnotationTable.from_hail(
            hl.read_table(paths["table"]), provenance=result.provenance
        ),
        paths["tsv"],
    )
    result.dictionary.to_csv(
        paths["dictionary"], sep="\t", index=False, lineterminator="\n"
    )
    Path(paths["manifest"]).write_text(
        json.dumps(result.manifest, indent=2, default=str) + "\n"
    )
    try:
        hl.current_backend().fs.rmtree(result.work_path)
    except Exception as exc:  # noqa: BLE001 -- leftover scratch is not worth failing for
        logger.warning("could not remove %s: %s", result.work_path, exc)
    return paths


# --- helpers ---------------------------------------------------------------------------


def check_inputs(spec, inputs, gene_matrix, gene_axes, gene_lists) -> None:
    """Every disagreement between spec, inputs and gene matrix, before any Hail job.

    Reads only table schemas. ``build_training_set`` runs it first; the CLI runs it on
    its own before loading anything, to report a mistake as a usage error.

    Raises
    ------
    ValueError
        Naming the disagreement.
    KeyError
        When a field path the spec or the builder relies on is not in its table.
    """
    if spec.label is None or spec.gene is None:
        raise ValueError(
            f"spec {spec.name!r} needs a label: and a gene: section to build a training set"
        )
    if spec.label.source not in LABEL_SOURCES:
        raise ValueError(
            f"label source {spec.label.source!r} is not supported; "
            f"supported: {', '.join(LABEL_SOURCES)}"
        )
    _same_names("inputs", spec.input_names(), inputs)
    _same_names("gene lists", [g.column for g in spec.gene_lists], gene_lists)
    for column, symbols in gene_lists.items():
        if not list(symbols):
            raise ValueError(f"gene list {column!r} is empty")
    for field in ("entrez_id", "ensembl_gene_id"):
        field_expr(inputs[spec.gene.hgnc_input], field)
    _check_gene_matrix(spec, gene_matrix, gene_axes)
    _check_variant_sources(spec, inputs)


def _check_gene_matrix(spec, gene_matrix, gene_axes) -> None:
    if list(gene_matrix.key) != ["gene_id"]:
        raise ValueError(
            f"the gene matrix must be keyed by gene_id (an annotate compose output); "
            f"it is keyed by {list(gene_matrix.key)}"
        )
    matrix_fields = set(gene_matrix.row)
    wanted = [c.name for c in spec.gene.spine_columns]
    for entry in spec.layer1:
        if entry.axis not in gene_axes:
            raise ValueError(
                f"axis {entry.axis!r} is not in the gene matrix's manifest; was the "
                "matrix composed from this spec?"
            )
        wanted += [*gene_axes[entry.axis], f"{entry.axis}_present"]
    missing = [c for c in wanted if c not in matrix_fields]
    if missing:
        raise ValueError(f"gene matrix lacks column(s): {', '.join(missing)}")
    taken = set(METADATA_COLUMNS) | {g.column for g in spec.gene_lists}
    taken |= {c.name for v in spec.variant for c in v.columns}
    clash = sorted(set(wanted) & taken)
    if clash:
        raise ValueError(
            "gene matrix column(s) collide with metadata, gene-list or variant "
            f"columns: {', '.join(clash)}"
        )


def _check_variant_sources(spec, inputs) -> None:
    for v in spec.variant:
        src = inputs[v.input]
        if v.join == "allele":
            if list(src.key) != ["locus", "alleles"]:
                raise ValueError(
                    f"variant source {v.input!r} joins on alleles but is keyed by "
                    f"{list(src.key)}, not [locus, alleles]"
                )
            for c in v.columns:
                for path in c.fields or (c.field,):
                    field_expr(src, path)
            if v.gene_field is not None:
                field_expr(src, v.gene_field)
        elif v.join == "interval":
            key = list(src.key)
            if len(key) != 1 or not isinstance(src[key[0]].dtype, hl.tinterval):
                raise ValueError(
                    f"variant source {v.input!r} joins on intervals but is keyed by "
                    f"{key}, not one interval field"
                )
        else:  # codon: the fields codon_sites reads
            for field in ("chrom", "codon_start", "codon_end"):
                field_expr(src, field)


def _same_names(what: str, expected: Sequence[str], given: Mapping) -> None:
    missing = sorted(set(expected) - set(given))
    unknown = sorted(set(given) - set(expected))
    if missing or unknown:
        parts = []
        if missing:
            parts.append("missing: " + ", ".join(missing))
        if unknown:
            parts.append("not in the spec: " + ", ".join(unknown))
        raise ValueError(f"{what} do not match the spec ({'; '.join(parts)})")


def _plan_columns(spec, gene_axes, score_provenance, input_provenance) -> list[_Column]:
    """Every column in output order (after the key), with what its tier derives from."""
    columns = [_metadata_column(name, spec) for name in _ROW_METADATA]
    columns += [
        _Column(
            gl.column, "gene_list", "gene", "gene list", gl.description, feature=False
        )
        for gl in spec.gene_lists
    ]
    for v in spec.variant:
        declared_scores = score_provenance.get(v.source, {})
        for c in v.columns:
            fields = c.fields or ((c.field,) if c.field is not None else ())
            trained = _trained_on(declared_scores, fields)
            columns.append(
                _Column(
                    name=c.name,
                    origin="variant",
                    level="variant" if v.join == "allele" else "protein",
                    source=v.source,
                    description=c.description
                    or f"{v.source} {c.field or ', '.join(c.fields) or v.join}",
                    input=v.input,
                    trained_on=trained,
                    declared=c.tier,
                )
            )
    for sc in spec.gene.spine_columns:
        columns.append(
            _Column(
                sc.name,
                "matrix",
                "gene",
                "gene spine",
                sc.description or sc.name,
                declared=sc.tier,
            )
        )
    for entry in spec.layer1:
        if entry.matrix is not None:
            what = f"per-{entry.matrix.group_axis} {'/'.join(entry.matrix.stats) or 'value'}"
        else:
            what = "value"
        # An axis aggregating per-variant scores (key: variant) carries those scores'
        # training provenance into its per-gene columns, so a declared tier cannot hide
        # circularity there either; its row-count column has none and keeps the axis tier.
        score_fields = {}
        if entry.aggregate is not None:
            score_fields = {
                output_name(sc.name, stat): sc.column
                for sc in entry.aggregate.scores
                for stat in sc.stats
            }
        declared_scores = score_provenance.get(entry.source, {})
        for col in gene_axes[entry.axis]:
            columns.append(
                _Column(
                    col,
                    "matrix",
                    "gene",
                    entry.source,
                    f"{entry.source} {what}: {col}",
                    trained_on=_trained_on(declared_scores, (score_fields[col],))
                    if col in score_fields
                    else None,
                    declared=entry.tier,
                )
            )
        columns.append(
            _Column(
                f"{entry.axis}_present",
                "matrix",
                "gene",
                entry.source,
                f"Gene has a row in axis {entry.axis} (false = absent from the source, "
                "not a measured zero)",
                declared=_PRESENT_TIER if entry.tier is not None else None,
            )
        )
    return columns


def _gene_list_report(
    symbols: Sequence[str], mapping: Mapping[str, str | None]
) -> dict:
    """How a gene list mapped onto the spine: counts that reconcile, every unmapped symbol.

    ``n_symbols`` distinct symbols = ``n_symbols_mapped`` + ``n_unmapped``; ``n_genes`` is
    the distinct gene ids they map to (two symbols can name one gene). The mapper's own
    report names at most 20 unmapped symbols, so the list is rebuilt from the mapping.
    """
    unmapped = sorted({s for s in symbols if mapping.get(s) is None})
    mapped = {s for s in symbols if mapping.get(s) is not None}
    return {
        "n_symbols": len(set(symbols)),
        "n_symbols_mapped": len(mapped),
        "n_genes": len({mapping[s] for s in mapped}),
        "n_unmapped": len(unmapped),
        "unmapped": unmapped,
    }


def _trained_on(declared_scores, fields) -> frozenset[str] | None:
    """Union of the declared ``trained_on`` of ``fields``; ``None`` if any is undeclared.

    A column built from several scores (``fields`` + ``reduce``) is as circular as its
    most circular input, and one undeclared input makes the whole column undeclared.
    """
    if not fields:
        return None
    sources: set[str] = set()
    for f in fields:
        declared = declared_scores.get(f)
        if declared is None:
            return None
        sources.update(declared)
    return frozenset(sources)


def _preferring_source(spec: FeatureSpec):
    """The variant source whose gene_field breaks multi-gene ties, if any."""
    return next((v for v in spec.variant if v.gene_field is not None), None)


def _metadata_column(name: str, spec: FeatureSpec) -> _Column:
    label = spec.label
    preferring = _preferring_source(spec)
    gene_rule = "the first GENEINFO gene"
    if preferring is not None:
        gene_rule = (
            f"the first GENEINFO gene {preferring.source} also assigns "
            f"({preferring.gene_field}), else the first GENEINFO gene"
        )
    described = {
        "locus": ("variant", "Variant locus"),
        "alleles": ("variant", "Reference and alternate allele"),
        "vid": ("variant", "chrom:pos:ref:alt"),
        "gene_id": (
            "gene",
            f"Ensembl gene id (GENEINFO Entrez id via HGNC) of {gene_rule}",
        ),
        "gene_symbol": ("gene", "GENEINFO symbol of the chosen gene"),
        "clinvar_clnsig": ("variant", "CLNSIG as released"),
        "clinvar_review_status": ("variant", "CLNREVSTAT as released"),
        "clinvar_review_stars": ("variant", "Review stars from CLNREVSTAT (0-4)"),
        "consequence": ("variant", "MC consequence terms"),
        "label": (
            "variant",
            f"1 = {', '.join(label.pathogenic)}; 0 = {', '.join(label.benign)} "
            "(primary CLNSIG)",
        ),
    }
    level, description = described[name]
    return _Column(
        name,
        "metadata",
        level,
        label.source,
        description,
        feature=False,
        input=label.input,
    )


def _source_version(input_provenance, input_name) -> str:
    prov = input_provenance.get(input_name) if input_name else None
    return "" if prov is None else str(prov.source_fingerprint)


def _str_keys(counter) -> dict:
    return {str(k): v for k, v in sorted(counter.items(), key=lambda kv: str(kv[0]))}


def _jsonable(value):
    return json.loads(json.dumps(value, default=str))


def _hvantk_version() -> str:
    try:
        return importlib.metadata.version("hvantk")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"
