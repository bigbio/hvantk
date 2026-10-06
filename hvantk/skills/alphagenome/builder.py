"""Hail Table builder for AlphaGenome variant-effect scores.

hvantk never calls AlphaGenome. It ingests the scores a user has already
produced with the AlphaGenome SDK (through the API or local weights):
``score_variant(...)`` -> ``variant_scorers.tidy_scores(...)`` -> parquet. See
SKILL.md s 2 for the producer snippet.

Input contract
--------------
``parsed_input`` is a directory holding one or more ``*.parquet`` files, or a
single parquet file. Only ``*.parquet`` files directly inside the directory are
read, so a NOTICE.md or a script kept beside them is ignored. Each file is the
SDK's ``tidy_scores()`` long format: one row per variant x scorer x track (x gene
or junction, for gene-level scorers). Columns used:

- required in every file: ``variant_id``, ``output_type``, ``variant_scorer``,
  ``track_name``, ``gene_id`` (null for the track-level scorers), ``raw_score``,
  ``quantile_score``;
- ``ontology_curie``, required only by the ``ontology_curies`` filter.

Each file is read with its own schema and cast to the same types (strings,
float64 scores) before the files are combined, so files may differ in physical
types: a ``gene_id`` stored as parquet ``null`` (a file holding only track-level
scorers), or float32 next to float64 scores.

``variant_id`` is the SDK's ``chrom:pos:ref>alt`` string, 1-based, on GRCh38
``chr`` contigs (e.g. ``chr3:39408741:T>C``). All other columns are ignored.

Output
------
One row per variant, keyed by ``(locus, alleles)``, carrying ``variant_id`` and
one struct per known scorer (``shared.constants.SCORER_FIELDS``), so scorers are
never mixed. Each struct summarises that scorer's rows for the variant:

- ``top_raw``, ``top_quantile``, ``top_track``, ``top_gene_id``: the row with the
  largest |quantile_score|, sign kept. Ties go to the smallest ``track_name``,
  then ``gene_id``, ``quantile_score`` and ``raw_score``, so the pick is stable.
- ``max_abs_raw``: the largest |raw_score|.
- ``n_rows``: the number of rows aggregated.

A struct is missing when the variant has no row for that scorer. Rows whose
``raw_score`` or ``quantile_score`` is missing or NaN are dropped first.

The rows are checked with Spark before anything reaches Hail. ``ValueError`` is
raised for a missing required column, an unknown ``output_types`` value, a
``variant_scorer`` not in ``SCORER_FIELDS`` anywhere in the input (whatever the
filters), a ``variant_id`` not in ``chrom:pos:ref>alt`` form, and filters that
leave no row to aggregate.
"""

from __future__ import annotations

from pathlib import Path

import hail as hl

from hvantk.skills.alphagenome.shared.constants import OUTPUT_TYPES, SCORER_FIELDS

SCHEMA_ID = "alphagenome-v2"
_STRING_COLUMNS = (
    "variant_id",
    "output_type",
    "variant_scorer",
    "track_name",
    "gene_id",
)
_SCORE_COLUMNS = ("raw_score", "quantile_score")
REQUIRED_COLUMNS = _STRING_COLUMNS + _SCORE_COLUMNS
#: ``str(genome.Variant)``, the SDK's default variant format: ``chrom:pos:ref>alt``
#: with bases from ACGTN (an allele may be empty).
VARIANT_ID_PATTERN = r"^[^:]+:[0-9]+:[ACGTN]*>[ACGTN]*$"


def _parquet_paths(input_path) -> list[str]:
    """Absolute paths of the input parquet files (Spark wants absolute paths)."""
    path = Path(input_path)
    files = [path] if path.is_file() else sorted(path.glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"No *.parquet file found at {input_path}")
    return [str(f.resolve()) for f in files]


def _as_list(value) -> list[str]:
    """``--plugin-arg`` passes a single value as a plain string, not a list."""
    return [value] if isinstance(value, str) else list(value)


def _read_scores(paths: list[str], need_ontology: bool):
    """Read the parquet files into one Spark DataFrame of the columns used.

    Spark reading several files at once applies one file's schema to all of them:
    a column typed differently in another file (``null`` vs string, float32 vs
    float64) then fails at read time, and a column the schema-giving file lacks
    reads as null everywhere. So each file is read on its own, checked, and cast
    to the same types, and the results are combined with ``unionByName``. The
    selection also keeps the ~18 unused columns, among them ``Assay title`` (with
    a space), out of Hail.
    """
    from functools import reduce

    from pyspark.sql import SparkSession
    from pyspark.sql import functions as F

    from hvantk.core.utils.hail_context import init_hail

    init_hail()  # Spark must be Hail's session, so Hail starts first
    spark = SparkSession.builder.getOrCreate()

    strings = _STRING_COLUMNS + (("ontology_curie",) if need_ontology else ())
    frames = []
    for path in paths:
        sdf = spark.read.parquet(path)
        missing = sorted(set(strings + _SCORE_COLUMNS) - set(sdf.columns))
        if missing:
            raise ValueError(
                f"{path} is missing required AlphaGenome column(s): "
                f"{', '.join(missing)}"
            )
        frames.append(
            sdf.select(
                *[F.col(c).cast("string").alias(c) for c in strings],
                *[F.col(c).cast("double").alias(c) for c in _SCORE_COLUMNS],
            )
        )
    return reduce(lambda a, b: a.unionByName(b), frames)


def _checked_rows(sdf, output_types, ontology_curies):
    """Check the rows on the Spark side, then keep the ones to aggregate.

    One Spark pass before anything reaches Hail, so a bad input fails before the
    shuffle. Each filter step is counted, so an empty result says which filter
    removed everything.
    """
    from pyspark.sql import functions as F

    def is_number(column):
        # NaN is a value, not a null, so it needs its own test.
        return F.col(column).isNotNull() & ~F.isnan(column)

    # Each filter as (label, condition to pass it and every filter before it).
    steps = [
        ("with numeric scores", is_number("raw_score") & is_number("quantile_score"))
    ]
    if output_types is not None:
        keep = steps[-1][1] & F.col("output_type").isin(output_types)
        steps.append((f"with output_types {output_types}", keep))
    if ontology_curies is not None:
        keep = steps[-1][1] & F.col("ontology_curie").isin(ontology_curies)
        steps.append((f"with ontology_curies {ontology_curies}", keep))

    variant_id = F.col("variant_id")
    bad_id = variant_id.isNull() | ~variant_id.rlike(VARIANT_ID_PATTERN)
    stats = sdf.agg(
        F.collect_set("variant_scorer").alias("scorers"),
        F.first(
            F.when(bad_id, F.coalesce(variant_id, F.lit("<null>"))), ignorenulls=True
        ).alias("bad_id"),
        F.count(F.lit(1)).alias("rows"),
        *[
            F.count(F.when(kept, True)).alias(f"step{i}")
            for i, (_, kept) in enumerate(steps)
        ],
    ).first()

    unknown = sorted(set(stats["scorers"]) - set(SCORER_FIELDS))
    if unknown:
        raise ValueError(
            "Unknown AlphaGenome variant_scorer value(s): "
            + "; ".join(unknown)
            + ". Add each to SCORER_FIELDS in "
            "hvantk/skills/alphagenome/shared/constants.py once its output is "
            "understood (see SKILL.md s 8)."
        )
    if stats["bad_id"] is not None:
        raise ValueError(
            f"variant_id {stats['bad_id']!r} is not in the SDK's chrom:pos:ref>alt "
            "form (e.g. chr3:39408741:T>C)"
        )
    counts = [("rows", stats["rows"])]
    counts += [(label, stats[f"step{i}"]) for i, (label, _) in enumerate(steps)]
    if counts[-1][1] == 0:
        raise ValueError(
            "No AlphaGenome rows left to aggregate: "
            + ", ".join(f"{n:,} {label}" for label, n in counts)
        )
    return sdf.filter(steps[-1][1])


def build_alphagenome_predictions(
    parsed_input,
    ctx,
    *,
    reference_genome: str = "GRCh38",
    output_types=None,
    ontology_curies=None,
):
    """Plugin builder: AlphaGenome tidy-scores parquet -> AnnotationTable.

    Parameters
    ----------
    parsed_input : str | Path
        Directory of ``tidy_scores()`` parquet files, or one parquet file.
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    reference_genome : str
        Reference genome of the ``locus`` key. AlphaGenome scores human
        variants on GRCh38.
    output_types : list[str] | str | None
        Keep only rows of these ``output_type`` values (e.g. ``["RNA_SEQ"]``);
        the other fields come out missing. ``None`` keeps all.
    ontology_curies : list[str] | str | None
        Keep only rows of tracks with these ``ontology_curie`` values before
        aggregating, e.g. heart: ``["UBERON:0006566", "UBERON:0006631"]``.
        ``None`` keeps all.
    """
    from hvantk.core.models import AnnotationTable

    if output_types is not None:
        output_types = _as_list(output_types)
        unknown = sorted(set(output_types) - OUTPUT_TYPES)
        if unknown:
            raise ValueError(
                f"Unknown output_types {unknown}; known: {sorted(OUTPUT_TYPES)}"
            )
    if ontology_curies is not None:
        ontology_curies = _as_list(ontology_curies)

    sdf = _read_scores(
        _parquet_paths(parsed_input), need_ontology=ontology_curies is not None
    )
    ht = hl.Table.from_spark(_checked_rows(sdf, output_types, ontology_curies))

    # One summary per (variant, scorer).
    top = hl.agg.take(
        hl.struct(
            top_raw=ht.raw_score,
            top_quantile=ht.quantile_score,
            top_track=ht.track_name,
            top_gene_id=ht.gene_id,
        ),
        1,
        ordering=hl.tuple(
            [
                -hl.abs(ht.quantile_score),
                ht.track_name,
                ht.gene_id,
                ht.quantile_score,
                ht.raw_score,
            ]
        ),
    )[0]
    pairs = ht.group_by(ht.variant_id, ht.variant_scorer).aggregate(
        summary=top.annotate(
            max_abs_raw=hl.agg.max(hl.abs(ht.raw_score)),
            n_rows=hl.agg.count(),
        )
    )

    # Pivot to one row per variant; a scorer the variant lacks stays missing.
    field = hl.literal(SCORER_FIELDS)[pairs.variant_scorer]
    ht = pairs.group_by(pairs.variant_id).aggregate(
        summaries=hl.dict(hl.agg.collect(hl.tuple([field, pairs.summary])))
    )
    parts = ht.variant_id.split(":")
    ht = ht.key_by(
        locus=hl.locus(parts[0], hl.int32(parts[1]), reference_genome=reference_genome),
        alleles=parts[2].split(">"),
    )
    ht = ht.select(
        "variant_id", **{f: ht.summaries.get(f) for f in SCORER_FIELDS.values()}
    )

    return AnnotationTable.from_hail(ht, provenance=ctx.provenance(schema_id=SCHEMA_ID))
