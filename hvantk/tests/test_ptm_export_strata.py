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


def test_ptm_population_excludes_missing_af_from_the_taken_arrays(
    tmp_path, hail_session
):
    """ptm_population's per-category AF sample used ``hl.agg.take(af, 5000)``, which
    keeps missing values (sorted last) rather than dropping them; ``float(x)`` on that
    array then raised TypeError the first time a PTM-site variant had no AF.
    """
    import hail as hl

    from hvantk.algorithms.ptm.analysis import ptm_population

    ptm_ht = hl.Table.parallelize(
        [
            {
                "locus": hl.Locus("chr1", 100, "GRCh38"),
                "codon_start": 100,
                "codon_end": 102,
                "ptm_category": "phosphorylation",
            },
        ],
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            codon_start=hl.tint32,
            codon_end=hl.tint32,
            ptm_category=hl.tstr,
        ),
        key="locus",
    )

    gnomad_ht = hl.Table.parallelize(
        [
            # PTM site, AF present
            {"locus": hl.Locus("chr1", 100, "GRCh38"), "AF": 0.01},
            # PTM site, AF missing -- must not crash the float() conversion
            {"locus": hl.Locus("chr1", 101, "GRCh38"), "AF": None},
            # non-PTM (far outside the flanking window)
            {"locus": hl.Locus("chr1", 500, "GRCh38"), "AF": 0.05},
        ],
        hl.tstruct(locus=hl.tlocus("GRCh38"), AF=hl.tfloat64),
        key="locus",
    )

    result = ptm_population(gnomad_ht, ptm_ht, str(tmp_path))

    assert result.n_ptm_site == 2
    assert result.ptm_site_afs == [0.01]
