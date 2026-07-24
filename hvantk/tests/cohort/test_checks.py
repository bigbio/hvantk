"""Cohort manifest validation against its on-disk table.

Moved out of the CLI (``hvantk.tools.cohort.cohort_cli``) into
``hvantk.algorithms.cohort.checks`` so a Python caller of ``attach()`` can run the same
checks the CLI runs, without importing Click. Pure Python -- no Hail -- so these run in
the default fast suite.
"""
import gzip

import pytest

from hvantk.algorithms.cohort.checks import (
    check_declared_columns_exist,
    check_key_column_exists,
    detect_delimiter,
    read_header,
)
from hvantk.algorithms.cohort.spec import CohortManifest, CohortPrior


def _manifest(key="gene_id", key_column=None, prior_col="minp", axes=()):
    return CohortManifest(
        name="demo",
        key=key,
        key_column=key_column,
        table="/unused.tsv",
        prior=CohortPrior(column=prior_col, direction="lower_is_better"),
        cohort_axes=axes,
    )


def test_read_header_reads_a_comma_delimited_table(tmp_path):
    p = tmp_path / "cohort.csv"
    p.write_text("gene_id,minp,n_case_var\nENSG1,0.01,5\n")
    assert read_header(str(p)) == ["gene_id", "minp", "n_case_var"]


def test_read_header_reads_a_tab_delimited_table_named_tsv(tmp_path):
    p = tmp_path / "cohort.tsv"
    p.write_text("gene_id\tminp\tn_case_var\nENSG1\t0.01\t5\n")
    assert read_header(str(p)) == ["gene_id", "minp", "n_case_var"]


def test_read_header_sniffs_tab_delimiter_from_an_unusual_extension(tmp_path):
    """finding 3: the old suffix whitelist (`.tsv`/`.txt`) mis-delimited any tab file
    with any other extension -- `.tab`, or no extension at all -- producing a
    misleading "table has no column 'gene'" error even though the column is present,
    just comma-split into one field with embedded tabs. Sniffing the header line for a
    literal tab is extension-agnostic."""
    p = tmp_path / "cohort.tab"
    p.write_text("gene_id\tminp\tn_case_var\nENSG1\t0.01\t5\n")
    assert read_header(str(p)) == ["gene_id", "minp", "n_case_var"]


def test_read_header_handles_a_gzipped_tsv_table(tmp_path):
    """finding 3: `Path("x.tsv.gz").suffix == ".gz"`, so the delimiter picked ','
    (wrong) and the file was opened in text mode over raw gzip bytes -- a bare
    UnicodeDecodeError on the gzip magic byte 0x8b, not a clean validation error."""
    p = tmp_path / "cohort.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("gene_id\tminp\tn_case_var\nENSG1\t0.01\t5\n")
    assert read_header(str(p)) == ["gene_id", "minp", "n_case_var"]


def test_read_header_raises_for_a_missing_table(tmp_path):
    with pytest.raises(ValueError, match="not found"):
        read_header(str(tmp_path / "nope.tsv"))


def test_read_header_raises_for_an_empty_table(tmp_path):
    p = tmp_path / "cohort.tsv"
    p.write_text("")
    with pytest.raises(ValueError, match="empty"):
        read_header(str(p))


def test_detect_delimiter_agrees_with_read_header_for_gzip(tmp_path):
    """attach_cmd's hl.import_table and validate's read_header must never disagree --
    this is the same sniff, exposed separately because import_table needs the raw
    delimiter character rather than the parsed header."""
    p = tmp_path / "cohort.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("gene_id\tminp\n")
    assert detect_delimiter(str(p)) == "\t"


def test_detect_delimiter_for_a_comma_table(tmp_path):
    p = tmp_path / "cohort.csv"
    p.write_text("gene_id,minp\n")
    assert detect_delimiter(str(p)) == ","


def test_check_key_column_exists_raises_when_missing():
    m = _manifest(key="symbol", key_column="gene")
    with pytest.raises(ValueError, match="gene"):
        check_key_column_exists(m, ["symbol", "minp"])


def test_check_key_column_exists_passes_when_present():
    m = _manifest(key="symbol", key_column="gene")
    check_key_column_exists(m, ["gene", "minp"])


def test_check_declared_columns_exist_raises_when_missing():
    m = _manifest(axes=())
    with pytest.raises(ValueError, match="minp"):
        check_declared_columns_exist(m, ["gene_id"])


def test_check_declared_columns_exist_passes_when_present():
    m = _manifest(axes=())
    check_declared_columns_exist(m, ["gene_id", "minp"])
