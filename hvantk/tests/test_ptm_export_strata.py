"""#307: export_ptm_strata makes two passes (one export per stratum), not three."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.hail


def test_export_ptm_strata_counts_come_from_the_exported_files(
    tmp_path, hail_session, monkeypatch
):
    import hail as hl

    from hvantk.algorithms.ptm import analysis

    ht = hl.Table.parallelize(
        [
            {
                "locus": hl.Locus("chr1", 100, "GRCh38"),
                "alleles": ["A", "T"],
                "is_ptm_site": True,
                "is_ptm_proximal": False,
            },
            {
                "locus": hl.Locus("chr1", 200, "GRCh38"),
                "alleles": ["G", "C"],
                "is_ptm_site": False,
                "is_ptm_proximal": True,
            },
            {
                "locus": hl.Locus("chr2", 300, "GRCh38"),
                "alleles": ["T", "A"],
                "is_ptm_site": False,
                "is_ptm_proximal": False,
            },
        ],
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            is_ptm_site=hl.tbool,
            is_ptm_proximal=hl.tbool,
        ),
        key=["locus", "alleles"],
    )
    # The aggregate is the third pass the issue is about: it must be gone.
    monkeypatch.setattr(
        hl.Table,
        "aggregate",
        lambda *a, **k: pytest.fail("export_ptm_strata still aggregates"),
    )

    paths = analysis.export_ptm_strata(ht, str(tmp_path))

    assert set(paths) == {"ptm", "non_ptm"}
    ptm_lines = (tmp_path / "ptm_variants.txt").read_text().splitlines()
    non_ptm_lines = (tmp_path / "non_ptm_variants.txt").read_text().splitlines()
    assert sorted(ptm_lines) == ["chr1:100:A:T", "chr1:200:G:C"]
    assert non_ptm_lines == ["chr2:300:T:A"]
