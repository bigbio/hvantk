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
  ``track_name``, ``track_strand``, ``gene_id`` (null for the track-level
  scorers), ``raw_score``, ``quantile_score``;
- ``junction_Start`` and ``junction_End``, read when present (null except on
  SPLICE_JUNCTIONS rows);
- ``ontology_curie``, required only by the ``ontology_curies`` filter.

Each file is read with its own schema and cast to the same types (strings,
float64 scores, int64 junction coordinates) before the files are combined, so
files may differ in physical types: a ``gene_id`` stored as parquet ``null`` (a
file holding only track-level scorers), or float32 next to float64 scores.

``variant_id`` is the SDK's ``chrom:pos:ref>alt`` string, 1-based, on GRCh38
``chr`` contigs (e.g. ``chr3:39408741:T>C``). All other columns are ignored.

A tidy-scores row is identified by ``ROW_KEY``: variant, scorer, track, strand,
gene and junction. The same key twice means the same scores were passed twice
(a file given twice, or two scoring runs mixed) and would be double-counted, so
the build fails. On a real 16.6M-row run every row has its own key, while
leaving ``track_strand``, the junction columns or ``gene_id`` out of the key
makes hundreds of thousands to millions of distinct rows collide.

Output
------
One row per variant that still has rows after the filters, keyed by
``(locus, alleles)``, carrying ``variant_id`` and one struct per known scorer
(``shared.constants.SCORER_FIELDS``), so scorers are never mixed. Each struct
summarises that scorer's rows for the variant:

- ``top_raw``, ``top_quantile``, ``top_track``, ``top_gene_id``: the row with the
  largest |quantile_score|, sign kept, among the rows that have a quantile score
  (missing when none has one). Ties go to the smallest ``track_name``, then
  ``gene_id``, ``quantile_score`` and ``raw_score``, so the pick is stable.
- ``max_abs_raw``: the largest |raw_score| (missing when no row has one).
- ``n_rows``: the number of rows aggregated.

A struct is missing when the variant has no row for that scorer. NaN scores are
treated as missing. A row with neither score is dropped; a row with only one of
the two still counts in ``n_rows`` and feeds the statistic it has a score for.
Dropped rows are reported in a warning with their count.

The rows are checked with Spark before anything reaches Hail. ``ValueError`` is
raised for a missing required column, an unknown ``output_types`` value, a
missing ``variant_scorer`` or one not in ``SCORER_FIELDS`` anywhere in the input
(whatever the filters), a ``variant_id`` not in ``chrom:pos:ref>alt`` form, a
``ROW_KEY`` that occurs more than once, and filters that leave no row to
aggregate.
"""

from __future__ import annotations

import logging
from pathlib import Path
from urllib.parse import urlsplit

import hail as hl

from hvantk.skills.alphagenome.shared.constants import OUTPUT_TYPES, SCORER_FIELDS

logger = logging.getLogger(__name__)

SCHEMA_ID = "alphagenome-v2"
_STRING_COLUMNS = (
    "variant_id",
    "output_type",
    "variant_scorer",
    "track_name",
    "track_strand",
    "gene_id",
)
_SCORE_COLUMNS = ("raw_score", "quantile_score")
REQUIRED_COLUMNS = _STRING_COLUMNS + _SCORE_COLUMNS
#: Read when present; null except on SPLICE_JUNCTIONS rows.
_JUNCTION_COLUMNS = ("junction_Start", "junction_End")
#: The identity of one tidy-scores row (see the module docstring).
ROW_KEY = (
    "variant_id",
    "variant_scorer",
    "track_name",
    "track_strand",
    "gene_id",
) + _JUNCTION_COLUMNS
#: ``str(genome.Variant)``, the SDK's default variant format: ``chrom:pos:ref>alt``
#: with non-empty alleles from ACGTN.
VARIANT_ID_PATTERN = r"^[^:]+:[0-9]+:[ACGTN]+>[ACGTN]+$"


def _parquet_paths(input_path, spark=None) -> list[str]:
    """Absolute paths of the input parquet files (Spark wants absolute paths)."""
    input_value = str(input_path)
    scheme = urlsplit(input_value).scheme
    if scheme in {"gs", "hdfs"}:
        if input_value.endswith(".parquet"):
            return [input_value]
        if spark is None:
            from pyspark.sql import SparkSession

            from hvantk.core.utils.hail_context import init_hail

            init_hail()
            spark = SparkSession.builder.getOrCreate()
        jvm = spark.sparkContext._jvm
        conf = spark.sparkContext._jsc.hadoopConfiguration()
        glob_path = jvm.org.apache.hadoop.fs.Path(
            input_value.rstrip("/") + "/*.parquet"
        )
        statuses = glob_path.getFileSystem(conf).globStatus(glob_path)
        if not statuses:
            raise FileNotFoundError(f"No *.parquet file found at {input_value}")
        return sorted(str(status.getPath().toString()) for status in statuses)

    path = Path(input_path)
    files = [path] if path.is_file() else sorted(path.glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"No *.parquet file found at {input_path}")
    return [str(f.resolve()) for f in files]


def _as_list(value) -> list[str]:
    """``--plugin-arg`` passes a single value as a plain string, not a list."""
    return [value] if isinstance(value, str) else list(value)


def _read_scores(input_path, need_ontology: bool):
    """Read the parquet files into one Spark DataFrame of the columns used.

    Spark reading several files at once applies one file's schema to all of them:
    a column typed differently in another file (``null`` vs string, float32 vs
    float64) then fails at read time, and a column the schema-giving file lacks
    reads as null everywhere. So each file is read on its own, checked, and cast
    to the same types, and the results are combined with ``unionByName``. The
    selection also keeps the unused columns, among them ``Assay title`` (with a
    space), out of Hail.
    """
    from functools import reduce

    from pyspark.sql import SparkSession
    from pyspark.sql import functions as F

    from hvantk.core.utils.hail_context import init_hail

    if urlsplit(str(input_path)).scheme not in {"gs", "hdfs"}:
        paths = _parquet_paths(input_path)
    else:
        paths = None
    init_hail()  # Spark must be Hail's session, so Hail starts first
    spark = SparkSession.builder.getOrCreate()
    if paths is None:
        paths = _parquet_paths(input_path, spark=spark)

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
                *[
                    (F.col(c) if c in sdf.columns else F.lit(None))
                    .cast("long")
                    .alias(c)
                    for c in _JUNCTION_COLUMNS
                ],
                *[F.col(c).cast("double").alias(c) for c in _SCORE_COLUMNS],
            )
        )
    return reduce(lambda a, b: a.unionByName(b), frames)


def _check_unique_rows(sdf) -> None:
    """Fail when a ``ROW_KEY`` occurs more than once (see the module docstring)."""
    from pyspark.sql import functions as F

    key_hash = F.sha2(F.to_json(F.struct(*[F.col(column) for column in ROW_KEY])), 256)
    hashed = sdf.withColumn("_row_key_hash", key_hash)
    possible_hashes = [
        row["_row_key_hash"]
        for row in (
            hashed.groupBy("_row_key_hash")
            .count()
            .filter(F.col("count") > 1)
            .select("_row_key_hash")
            .collect()
        )
    ]
    if not possible_hashes:
        return

    duplicated = (
        hashed.filter(F.col("_row_key_hash").isin(*possible_hashes))
        .groupBy(*ROW_KEY)
        .count()
        .filter(F.col("count") > 1)
    )
    example = duplicated.limit(1).collect()
    if example:
        row = example[0]
        raise ValueError(
            f"{duplicated.count():,} AlphaGenome score row(s) occur more than once "
            f"with the same {', '.join(ROW_KEY)}; for example "
            f"{row['variant_id']} / {row['variant_scorer']} / track "
            f"{row['track_name']!r} occurs {row['count']} times. The input holds "
            "the same scores twice: a file given twice, two scoring runs mixed, or "
            "the same variant scored twice in one run. The summaries would count "
            "them twice, so keep one copy of each. A parquet without "
            "junction_Start/junction_End also collides on its SPLICE_JUNCTIONS rows."
        )


def _checked_rows(sdf, output_types, ontology_curies):
    """Check the rows on the Spark side, then keep the ones to aggregate.

    Two Spark passes before anything reaches Hail, so a bad input fails before
    the shuffle: one validates the rows and counts each filter step, the other
    (``_check_unique_rows``) looks for a ``ROW_KEY`` that occurs twice. Rows
    dropped for lacking both scores are reported in a warning; filters that
    leave nothing raise, naming the step that removed everything.
    """
    from pyspark.sql import functions as F

    # NaN is a value, not a null: make it null so Hail's aggregators skip it.
    sdf = sdf.select(
        *[
            F.when(F.isnan(c), F.lit(None)).otherwise(F.col(c)).alias(c)
            if c in _SCORE_COLUMNS
            else F.col(c)
            for c in sdf.columns
        ]
    )
    has_raw = F.col("raw_score").isNotNull()
    has_quantile = F.col("quantile_score").isNotNull()

    # Each filter as (label, condition to pass it and every filter before it).
    steps = [("with a raw or quantile score", has_raw | has_quantile)]
    if output_types is not None:
        keep = steps[-1][1] & F.col("output_type").isin(output_types)
        steps.append((f"with output_types {output_types}", keep))
    if ontology_curies is not None:
        keep = steps[-1][1] & F.col("ontology_curie").isin(ontology_curies)
        steps.append((f"with ontology_curies {ontology_curies}", keep))
    kept = steps[-1][1]

    variant_id = F.col("variant_id")
    bad_id = variant_id.isNull() | ~variant_id.rlike(VARIANT_ID_PATTERN)
    stats = sdf.agg(
        F.collect_set("variant_scorer").alias("scorers"),
        F.count(F.when(F.col("variant_scorer").isNull(), True)).alias("no_scorer"),
        F.first(
            F.when(bad_id, F.coalesce(variant_id, F.lit("<null>"))), ignorenulls=True
        ).alias("bad_id"),
        F.count(F.lit(1)).alias("rows"),
        *[
            F.count(F.when(passed, True)).alias(f"step{i}")
            for i, (_, passed) in enumerate(steps)
        ],
        F.count(F.when(kept & ~has_quantile, True)).alias("raw_only"),
        F.count(F.when(kept & ~has_raw, True)).alias("quantile_only"),
    ).first()

    if stats["no_scorer"]:
        raise ValueError(
            f"{stats['no_scorer']:,} AlphaGenome row(s) have no variant_scorer"
        )
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

    _check_unique_rows(sdf)

    counts = [("rows", stats["rows"])]
    counts += [(label, stats[f"step{i}"]) for i, (label, _) in enumerate(steps)]
    if counts[-1][1] == 0:
        raise ValueError(
            "No AlphaGenome rows left to aggregate: "
            + ", ".join(f"{n:,} {label}" for label, n in counts)
        )
    no_score = stats["rows"] - stats["step0"]
    if no_score:
        logger.warning(
            "AlphaGenome: dropped %s of %s row(s) that have neither a raw nor a "
            "quantile score.",
            f"{no_score:,}",
            f"{stats['rows']:,}",
        )
    if stats["raw_only"] or stats["quantile_only"]:
        logger.warning(
            "AlphaGenome: %s aggregated row(s) have a raw score but no quantile "
            "score (they count in n_rows and max_abs_raw, not in the top_* pick); "
            "%s have a quantile score but no raw score (they count in n_rows and "
            "can be the top_* pick, with top_raw missing).",
            f"{stats['raw_only']:,}",
            f"{stats['quantile_only']:,}",
        )
    if len(steps) > 1:
        logger.info(
            "AlphaGenome filters: %s",
            ", ".join(f"{n:,} {label}" for label, n in counts),
        )
    return sdf.filter(kept).drop("track_strand", *_JUNCTION_COLUMNS)


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
        the other scorers' fields come out missing, and a variant with no row
        of these types is absent. ``None`` keeps all.
    ontology_curies : list[str] | str | None
        Keep only rows of tracks with these ``ontology_curie`` values before
        aggregating, e.g. heart: ``["UBERON:0006566", "UBERON:0006631"]``. A
        variant with no such track is absent. ``None`` keeps all.
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

    sdf = _read_scores(parsed_input, need_ontology=ontology_curies is not None)
    ht = hl.Table.from_spark(_checked_rows(sdf, output_types, ontology_curies))

    # One summary per (variant, scorer). Missing scores are skipped: the top_*
    # pick only considers rows with a quantile score, and max_abs_raw only rows
    # with a raw score.
    top = hl.agg.filter(
        hl.is_defined(ht.quantile_score),
        hl.agg.take(
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
        ),
    )
    pairs = ht.group_by(ht.variant_id, ht.variant_scorer).aggregate(
        top=top,
        max_abs_raw=hl.agg.max(hl.abs(ht.raw_score)),
        n_rows=hl.agg.count(),
    )
    pick = hl.or_missing(hl.len(pairs.top) > 0, pairs.top[0])
    pairs = pairs.select(
        summary=hl.struct(
            top_raw=pick.top_raw,
            top_quantile=pick.top_quantile,
            top_track=pick.top_track,
            top_gene_id=pick.top_gene_id,
            max_abs_raw=pairs.max_abs_raw,
            n_rows=pairs.n_rows,
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

    return AnnotationTable.from_hail(
        ht,
        provenance=ctx.provenance(
            schema_id=SCHEMA_ID,
            build_parameters={
                "output_types": output_types,
                "ontology_curies": ontology_curies,
            },
        ),
    )
