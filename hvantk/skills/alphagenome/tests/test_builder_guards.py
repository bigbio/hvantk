"""Guard tests for the AlphaGenome tidy-scores builder.

Each test feeds the builder a doctored copy of the real-subset fixture
(AlphaGenome Output Terms, non-commercial; see ``testdata/raw/alphagenome/NOTICE.md``)
and checks one documented rule: the checks that fail the build (unknown or missing
scorer, malformed id, duplicated rows, filters that leave nothing), missing-score
handling, input split across files of different float types, and the
``output_types`` / ``ontology_curies`` filters given as a plain string.
"""

from __future__ import annotations

import logging
from itertools import count
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from hvantk.skills.alphagenome.shared.constants import SCORER_FIELDS
from hvantk.tests._snapshot_utils import phase_b_snapshot_adapter

pytestmark = pytest.mark.hail

FIXTURE = next((Path(__file__).parent / "testdata/raw/alphagenome").glob("*.parquet"))
SCORER = {field: scorer for scorer, field in SCORER_FIELDS.items()}
CHR3, CHR6, CHRX = "chr3:39408741:T>C", "chr6:112216367:C>A", "chrX:153694448:T>G"
LOGGER = "hvantk.skills.alphagenome.builder"


def _base() -> pd.DataFrame:
    return pd.read_parquet(FIXTURE)


def _write(frame: pd.DataFrame, path: Path, float32: bool) -> None:
    """Write ``frame`` as parquet, keeping NaN scores as NaN (not null)."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    table = pa.Table.from_pandas(frame, preserve_index=False)
    for column in ("raw_score", "quantile_score"):
        values = frame[column].to_numpy(dtype="float32" if float32 else "float64")
        index = table.schema.get_field_index(column)
        # pa.array keeps NaN; Table.from_pandas would have turned it into null.
        table = table.set_column(index, column, pa.array(values))
    pq.write_table(table, path)


@pytest.fixture
def build(hail_session, tmp_path):
    """Build from one or more frames, one parquet file each; returns the table."""
    from hvantk.skills.alphagenome.builder import build_alphagenome_predictions

    builder = phase_b_snapshot_adapter(
        build_alphagenome_predictions, "alphagenome:predictions"
    )
    runs = count()

    def _build(*frames, float32=(), **params):
        run = next(runs)
        input_dir = tmp_path / f"input{run}"
        input_dir.mkdir()
        for i, frame in enumerate(frames):
            _write(frame, input_dir / f"part{i}.parquet", float32=i in float32)
        return builder(
            input_path=str(input_dir),
            output_path=str(tmp_path / f"output{run}.ht"),
            **params,
        )

    return _build


def _rows(ht) -> dict:
    return {row.variant_id: row for row in ht.collect()}


def _total_rows(rows: dict) -> int:
    """Rows aggregated over every variant and scorer."""
    return sum(
        getattr(row, field).n_rows
        for row in rows.values()
        for field in SCORER_FIELDS.values()
        if getattr(row, field) is not None
    )


def test_unknown_scorer_fails_even_when_filtered_out(build):
    df = _base()
    df.loc[df.index[df.output_type == "RNA_SEQ"][0], "variant_scorer"] = "MadeUp()"
    with pytest.raises(ValueError, match="Unknown AlphaGenome variant_scorer"):
        build(df, output_types=["SPLICE_SITES"])


def test_missing_scorer_fails(build):
    df = _base()
    df.loc[df.index[0], "variant_scorer"] = None
    with pytest.raises(ValueError, match="have no variant_scorer"):
        build(df)


def test_malformed_variant_id_fails(build):
    df = _base()
    df.loc[df.index[df.variant_id == CHR6][0], "variant_id"] = "chr6:112216367:C:A"
    with pytest.raises(ValueError, match="not in the SDK's chrom:pos:ref>alt form"):
        build(df)


@pytest.mark.parametrize("variant_id", ["chr3:39408741:>C", "chr3:39408741:T>"])
def test_empty_variant_allele_fails(build, variant_id):
    df = _base()
    df.loc[df.index[0], "variant_id"] = variant_id
    with pytest.raises(ValueError, match="not in the SDK's chrom:pos:ref>alt form"):
        build(df)


def test_parquet_paths_enumerate_cloud_files():
    from hvantk.skills.alphagenome.builder import _parquet_paths

    class HadoopPath:
        def __init__(self, value):
            self.value = value

        def getFileSystem(self, _conf):
            return filesystem

    class FileStatus:
        def __init__(self, path):
            self.path = path

        def getPath(self):
            return SimpleNamespace(toString=lambda: self.path)

    class FileSystem:
        def globStatus(self, path):
            assert path.value == "gs://bucket/scores/*.parquet"
            return [FileStatus("gs://bucket/scores/part-000.parquet")]

    filesystem = FileSystem()
    jvm = SimpleNamespace(
        org=SimpleNamespace(
            apache=SimpleNamespace(
                hadoop=SimpleNamespace(fs=SimpleNamespace(Path=HadoopPath))
            )
        )
    )
    spark = SimpleNamespace(
        sparkContext=SimpleNamespace(
            _jvm=jvm, _jsc=SimpleNamespace(hadoopConfiguration=lambda: object())
        )
    )

    assert _parquet_paths("gs://bucket/scores", spark=spark) == [
        "gs://bucket/scores/part-000.parquet"
    ]
    assert _parquet_paths("hdfs://namenode/data/part.parquet") == [
        "hdfs://namenode/data/part.parquet"
    ]


def test_same_scores_twice_fail(build):
    """A file given twice would double every n_rows."""
    df = _base()
    with pytest.raises(ValueError, match="occur more than once"):
        build(df, df)


def test_rows_differing_only_in_strand_or_junction_are_distinct(build):
    """``ROW_KEY`` includes the strand and junction columns, so these rows build."""
    from hvantk.skills.alphagenome.builder import ROW_KEY

    df = _base()
    key = list(ROW_KEY)
    existing = set(df[key].astype(str).itertuples(index=False, name=None))

    def copy_with_new_key(mask, column, change):
        for index in df.index[mask]:
            copy = df.loc[[index]].copy()
            copy[column] = change(copy[column].iloc[0])
            if tuple(copy[key].astype(str).iloc[0]) not in existing:
                return copy
        raise AssertionError(f"no row whose changed {column} gives a new key")

    other_strand = copy_with_new_key(
        df.output_type == "DNASE", "track_strand", lambda s: "-" if s == "+" else "+"
    )
    other_junction = copy_with_new_key(
        df.output_type == "SPLICE_JUNCTIONS", "junction_End", lambda end: end + 1
    )
    both = pd.concat([df, other_strand, other_junction])
    assert _total_rows(_rows(build(both))) == len(both)


def test_filters_that_leave_nothing_fail(build):
    with pytest.raises(ValueError, match="No AlphaGenome rows left to aggregate"):
        build(_base(), ontology_curies=["UBERON:9999999"])


def test_missing_scores_keep_what_they_have(build, caplog):
    caplog.set_level(logging.WARNING, logger=LOGGER)
    df = _base()
    splice = df.variant_scorer == SCORER["splice_sites"]
    # chr6: the donor row has both the largest |raw| and the largest |quantile|.
    donor = df.index[splice & (df.variant_id == CHR6) & (df.track_name == "donor")]
    df.loc[donor, "quantile_score"] = float("nan")
    # chrX: no SPLICE_SITES row has a quantile score.
    df.loc[splice & (df.variant_id == CHRX), "quantile_score"] = float("nan")
    # chr3: the donor row (the top |quantile|) has no raw score, and one row of
    # another scorer has neither score, so it is dropped.
    chr3_splice = splice & (df.variant_id == CHR3)
    df.loc[chr3_splice & (df.track_name == "donor"), "raw_score"] = float("nan")
    acceptor_raw = abs(df.loc[chr3_splice & (df.track_name == "acceptor"), "raw_score"])
    dropped = df.index[(df.variant_id == CHR3) & ~splice][0]
    df.loc[dropped, ["raw_score", "quantile_score"]] = float("nan")

    rows = _rows(build(df))

    chr3 = rows[CHR3].splice_sites
    assert chr3.n_rows == 2  # the donor row still counts
    assert chr3.top_track == "donor" and chr3.top_raw is None  # still the top pick
    assert chr3.max_abs_raw == pytest.approx(float(acceptor_raw.iloc[0]), rel=1e-6)
    chr6 = rows[CHR6].splice_sites
    assert chr6.max_abs_raw == pytest.approx(0.9982147, rel=1e-6)  # donor still counts
    assert chr6.n_rows == 2
    assert chr6.top_track == "acceptor"  # the only row with a quantile score
    chrx = rows[CHRX].splice_sites
    assert chrx.max_abs_raw == pytest.approx(0.9948883, rel=1e-6)
    assert chrx.n_rows == int((splice & (df.variant_id == CHRX)).sum())
    assert chrx.top_track is None and chrx.top_quantile is None
    assert _total_rows(rows) == len(df) - 1
    assert "dropped 1 of 150 row(s) that have neither" in caplog.text
    assert "have a raw score but no quantile score" in caplog.text
    assert "have a quantile score but no raw score" in caplog.text


def test_files_with_different_float_types_combine(build):
    df = _base()
    chr6, rest = df[df.variant_id == CHR6], df[df.variant_id != CHR6]
    rows = _rows(build(chr6, rest, float32=(1,)))
    assert set(rows) == {CHR3, CHR6, CHRX}
    assert rows[CHR6].splice_sites.max_abs_raw == pytest.approx(0.9982147, rel=1e-6)
    assert _total_rows(rows) == len(df)


def test_output_types_given_as_a_string(build):
    """``--plugin-arg output_types=SPLICE_SITES`` arrives as a plain string."""
    rows = _rows(build(_base(), output_types="SPLICE_SITES"))
    assert set(rows) == {CHR3, CHR6, CHRX}
    for row in rows.values():
        assert row.splice_sites is not None
        assert all(
            getattr(row, field) is None
            for field in SCORER_FIELDS.values()
            if field != "splice_sites"
        )


def test_ontology_curies_given_as_a_string(build):
    """``--plugin-arg ontology_curies=UBERON:0006566`` arrives as a plain string."""
    as_string = build(_base(), ontology_curies="UBERON:0006566").collect()
    as_list = build(_base(), ontology_curies=["UBERON:0006566"]).collect()
    assert as_string and as_string == as_list
