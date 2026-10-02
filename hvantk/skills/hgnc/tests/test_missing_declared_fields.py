"""Unit tests for the hgnc builder's missing-declared-field check.

Pure Python, no Hail needed. build_hgnc_gene_lookup's rename map used to
silently tolerate any HGNC_GENE_FIELDS key absent from the input header --
upstream dropped location_sortable and only the drift probe noticed (#381).
_missing_declared_fields() is the extracted computation that lets the
builder warn about this instead.
"""

from __future__ import annotations

from pathlib import Path

from hvantk.skills.hgnc.builder import _missing_declared_fields
from hvantk.skills.hgnc.shared.constants import HGNC_GENE_FIELDS

FIXTURE = Path("hvantk/skills/hgnc/tests/testdata/raw/hgnc/hgnc_test_sample.tsv")


def test_declared_fields_all_present_in_the_live_fixture_header():
    """HGNC_GENE_FIELDS must not declare a column the live dump no longer ships.

    hgnc_test_sample.tsv is derived from a real download, so its header is
    the closest thing to "what HGNC currently publishes" this suite can check
    without a network fetch. Catches the #381 case (location_sortable
    declared but dropped upstream) directly, instead of relying on the
    build-time warning or the drift probe noticing first.
    """
    header = FIXTURE.read_text().splitlines()[0].split("\t")
    missing = set(HGNC_GENE_FIELDS) - set(header)
    assert not missing, (
        f"declared in HGNC_GENE_FIELDS but absent from {FIXTURE}: {missing}"
    )


def test_missing_declared_fields_is_empty_when_header_has_everything():
    assert _missing_declared_fields(set(HGNC_GENE_FIELDS)) == []


def test_missing_declared_fields_reports_each_absent_key():
    row_fields = set(HGNC_GENE_FIELDS) - {"location", "omim_id"}
    assert set(_missing_declared_fields(row_fields)) == {"location", "omim_id"}


def test_missing_declared_fields_ignores_undeclared_extra_columns():
    """An extra header column not in HGNC_GENE_FIELDS is not "missing" -- it is
    simply unmapped, which the existing rename-map filtering already handles.
    """
    row_fields = set(HGNC_GENE_FIELDS) | {"some_new_upstream_column"}
    assert _missing_declared_fields(row_fields) == []
