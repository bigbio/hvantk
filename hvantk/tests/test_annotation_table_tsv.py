"""AnnotationTable exports as deterministic gzipped TSV (pandas tests need no Java)."""

from __future__ import annotations

import gzip
from datetime import datetime, timezone

import click
import pandas as pd
import pytest

from hvantk.core import io as core_io
from hvantk.core.io._errors import ArtifactTypeError
from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.provenance import Provenance


def _prov():
    return Provenance(
        plugin="t",
        dataset="t:rows",
        plugin_version="0.0",
        source_fingerprint="sha256:abc",
        schema_id="t-rows-v1",
        build_timestamp=datetime(2026, 5, 20, tzinfo=timezone.utc),
        builder_commit=None,
    )


def _table(df):
    return AnnotationTable.from_pandas(df, provenance=_prov())


def _text(path):
    with gzip.open(path, "rb") as fh:
        return fh.read().decode("utf-8")


def test_tsv_save_is_deterministic_and_joins_arrays(tmp_path):
    df = pd.DataFrame(
        {
            "gene": ["B", "A", "C", "D"],
            "tags": [["A", "B|C%"], [], None, {"y", "x"}],
            "flag": [True, False, True, False],
            "n": pd.array([5, None, 7, 0], dtype="Int64"),
            "score": [1.5, None, 2.0, 0.25],
            "note": ["é", "plain", None, ""],
        }
    )
    ann = _table(df)
    first, second = tmp_path / "a.tsv.gz", tmp_path / "b.tsv.gz"
    core_io.save(ann, first)
    core_io.save(ann, second)

    raw = first.read_bytes()
    assert raw == second.read_bytes()
    # gzip header: MTIME (bytes 4-8) zero and no FNAME flag, so no clock or path
    assert raw[4:8] == b"\0\0\0\0" and not raw[3] & 0x08
    assert first.with_name("a.tsv.gz.provenance.json").exists()
    # sorted by every field left to right; arrays "|"-joined (sets sorted) with
    # "%" and "|" percent-encoded inside an element, empty array and missing values
    # empty, bools true/false, ints plain decimal.
    assert _text(first) == (
        "gene\ttags\tflag\tn\tscore\tnote\n"
        "A\t\tfalse\t\t\tplain\n"
        "B\tA|B%7CC%25\ttrue\t5\t1.5\té\n"
        "C\t\ttrue\t7\t2.0\t\n"
        "D\tx|y\tfalse\t0\t0.25\t\n"
    )


@pytest.mark.parametrize("bad", ["a\tb", "a\nb", "a\rb"])
def test_tsv_save_refuses_a_tab_in_a_cell(tmp_path, bad):
    ann = _table(pd.DataFrame({"gene": ["G1", "G2"], "note": ["ok", bad]}))
    out = tmp_path / "t.tsv.gz"
    with pytest.raises(ValueError, match=r"note.*G2"):
        core_io.save(ann, out)
    assert not out.exists()


def test_tsv_cannot_be_loaded(tmp_path):
    out = tmp_path / "t.tsv.gz"
    core_io.save(_table(pd.DataFrame({"g": ["a"]})), out)
    with pytest.raises(ArtifactTypeError, match="export"):
        core_io.load(out)


def test_reprocess_accepts_tsv_gz_for_an_annotation_table():
    from types import SimpleNamespace

    from hvantk.tools.plugins.reprocess_cli import (
        _check_output_extension,
        _expected_extensions,
    )

    for backend, native in (("pandas", ".parquet"), ("hail", ".ht")):
        spec = SimpleNamespace(
            name="x:y",
            backend=backend,
            artifact_type=SimpleNamespace(__name__="AnnotationTable"),
        )
        assert _expected_extensions(spec) == (native, ".tsv.gz")
        _check_output_extension(spec, "out.tsv.gz")
        _check_output_extension(spec, "out" + native)
        with pytest.raises(click.UsageError):
            _check_output_extension(spec, "out.csv")
    vm = SimpleNamespace(
        name="x:y",
        backend="hail",
        artifact_type=SimpleNamespace(__name__="VariantMatrix"),
    )
    with pytest.raises(click.UsageError):
        _check_output_extension(vm, "out.tsv.gz")


def test_tsv_big_ints_in_an_all_numeric_table_keep_every_digit(tmp_path):
    big = 2**60 + 1
    ann = _table(pd.DataFrame({"n": [big + 2, big], "s": [0.5, 1.5]}))
    out = tmp_path / "t.tsv.gz"
    core_io.save(ann, out)
    assert _text(out) == f"n\ts\n{big}\t1.5\n{big + 2}\t0.5\n"


def test_tsv_nullable_int_column_prints_integers_and_an_empty_cell(tmp_path):
    ann = _table(
        pd.DataFrame({"g": ["a", "b", "c"], "n": pd.array([3, pd.NA, 12], "Int64")})
    )
    assert ann.schema["n"] == "int"
    out = tmp_path / "t.tsv.gz"
    core_io.save(ann, out)
    assert _text(out) == "g\tn\na\t3\nb\t\nc\t12\n"


def test_tsv_sorts_ties_by_the_next_field_with_missing_last(tmp_path):
    df = pd.DataFrame(
        {"a": ["x", "x", "x", "w"], "b": [None, "z", "y", "q"], "v": [1, 2, 3, 4]}
    )
    out = tmp_path / "t.tsv.gz"
    core_io.save(_table(df), out)
    assert _text(out) == "a\tb\tv\nw\tq\t4\nx\ty\t3\nx\tz\t2\nx\t\t1\n"


def test_tsv_float32_nan_is_a_missing_cell(tmp_path):
    import numpy as np

    df = pd.DataFrame({"g": ["a"], "x": np.array([np.nan], dtype="float32")})
    out = tmp_path / "t.tsv.gz"
    core_io.save(_table(df), out)
    assert _text(out) == "g\tx\na\t\n"


@pytest.mark.parametrize(
    "bad", [[["a"], "b"], {"name": "x"}, [{"name": "x", "url": "u"}]]
)
def test_tsv_refuses_a_nested_array_and_a_struct_or_dict(tmp_path, bad):
    ann = _table(pd.DataFrame({"gene": ["G1", "G2"], "tags": [["ok"], bad]}))
    out = tmp_path / "t.tsv.gz"
    with pytest.raises(ValueError, match=r"tags.*G2.*(array element|struct or dict)"):
        core_io.save(ann, out)
    assert not out.exists()
    assert not out.with_name(out.name + ".partial").exists()


def test_tsv_empty_array_missing_and_blank_element_write_the_same_cell(tmp_path):
    df = pd.DataFrame({"g": ["a", "b", "c"], "t": [[], None, [""]]})
    out = tmp_path / "t.tsv.gz"
    core_io.save(_table(df), out)
    assert _text(out) == "g\tt\na\t\nb\t\nc\t\n"


@pytest.mark.hail
def test_tsv_hail_table_sorts_loci_genomically_without_to_pandas(tmp_path, monkeypatch):
    import hail as hl

    def boom(*a, **k):
        raise AssertionError("the hail export must not collect via to_pandas")

    monkeypatch.setattr(hl.Table, "to_pandas", boom)
    schema = hl.tstruct(
        locus=hl.tlocus("GRCh38"),
        alleles=hl.tarray(hl.tstr),
        mc=hl.tarray(hl.tstr),
        af=hl.tfloat64,
        f32=hl.tfloat32,
        info=hl.tstruct(n=hl.tint64),
        lnk=hl.tarray(hl.tstruct(name=hl.tstr, url=hl.tstr)),
    )

    def row(contig, pos, alleles, mc, af, f32, n, lnk=()):
        locus = hl.Locus(contig, pos, "GRCh38")
        info = None if n is None else hl.Struct(n=n)
        lnk = [hl.Struct(name=a, url=b) for a, b in lnk]
        return dict(
            locus=locus, alleles=alleles, mc=mc, af=af, f32=f32, info=info, lnk=lnk
        )

    rows = [
        row(
            "chr10",
            5,
            ["A", "T"],
            ["SO:1|syn"],
            0.5,
            0.1,
            1,
            [("a", "u|1"), ("b", None)],
        ),
        row("chr1", 1000000, ["G", "C"], [], float("nan"), None, 2**60 + 1),
        row("chr1", 20, ["C", "G"], None, None, 1.5, None),
        row("chr2", 7, ["T", "A"], ["50%", "RORα"], 1e-5, 0.25, 3),
    ]
    ht = hl.Table.parallelize(
        rows, schema=schema, key=["locus", "alleles"], n_partitions=2
    )
    ann = AnnotationTable.from_hail(ht, provenance=_prov())
    first, second = tmp_path / "a.tsv.gz", tmp_path / "b.tsv.gz"
    core_io.save(ann, first)
    core_io.save(ann, second)
    assert first.read_bytes() == second.read_bytes()
    # genomic order (not chr1:1000000 < chr1:20 < chr10 < chr2); "|" and "%"
    # percent-encoded inside an element; non-ASCII text kept as UTF-8 (Hail's
    # text mode reads ISO-8859-1, which wrote "RORα" as "RORÎ±"); NaN empty;
    # float32 at float32 precision; the struct flattened to info.n with every digit
    # of an int64 kept; the array of structs split into lnk.name / lnk.url arrays
    # whose elements line up
    assert _text(first) == (
        "locus\talleles\tmc\taf\tf32\tinfo.n\tlnk.name\tlnk.url\n"
        "chr1:20\tC|G\t\t\t1.5\t\t\t\n"
        "chr1:1000000\tG|C\t\t\t\t1152921504606846977\t\t\n"
        "chr2:7\tT|A\t50%25|RORα\t1e-05\t0.25\t3\t\t\n"
        "chr10:5\tA|T\tSO:1%7Csyn\t0.5\t0.1\t1\ta|b\tu%7C1|\n"
    )
