"""Fast unit + CLI tests for the GWAS → effector coloc pipeline (gwas_coloc).

The end-to-end positive control (AF→MYOZ1 CONFIRMED) lives in the
network-marked test_gwas_coloc_myoz1_network.py (excluded from the fast run).
"""

import numpy as np
from click.testing import CliRunner

from hvantk.algorithms.qtlcascade.gwas_coloc import (
    coloc_abf_two_traits,
    run_locus_coloc,
    eqtl_catalogue_url,
)


# ---------------------------------------------------------------------------
# ABF kernel — correctness on synthetic data
# ---------------------------------------------------------------------------


def _signal(n, peak_idx, peak_z, se):
    z = np.linspace(-0.5, 0.5, n)  # mild deterministic noise
    z[peak_idx] = peak_z
    return z * se, np.full(n, se)


def test_coloc_abf_shared_signal_high_pp4():
    b1, s1 = _signal(100, 50, 8.0, 0.02)
    b2, s2 = _signal(100, 50, 8.0, 0.03)  # same peak variant -> shared
    res = coloc_abf_two_traits(b1, s1, b2, s2)
    assert res["H4"] > 0.9
    assert res["H3"] < 0.1
    assert res["lead_idx"] == 50


def test_coloc_abf_distinct_signal_high_pp3():
    b1, s1 = _signal(100, 50, 8.0, 0.02)
    b2, s2 = _signal(100, 10, 8.0, 0.03)  # different peak variant -> distinct
    res = coloc_abf_two_traits(b1, s1, b2, s2)
    assert res["H3"] > 0.9
    assert res["H4"] < 0.1


def test_coloc_abf_empty_returns_h0():
    res = coloc_abf_two_traits(np.array([]), np.array([]), np.array([]), np.array([]))
    assert res["H0"] == 1.0 and res["n_variants"] == 0


# ---------------------------------------------------------------------------
# Region driver — ranks the shared-signal gene above the distinct one
# ---------------------------------------------------------------------------


def _gwas_dict(n, peak_idx, peak_z, se=0.02):
    z = np.linspace(-0.5, 0.5, n)
    z[peak_idx] = peak_z
    return {(i + 1, "A", "G"): (z[i] * se, se, 1.0) for i in range(n)}


def _eqtl_recs(n, peak_idx, peak_z, se=0.03):
    z = np.linspace(-0.5, 0.5, n)
    z[peak_idx] = peak_z
    return [((i + 1, "A", "G"), z[i] * se, se, 1.0) for i in range(n)]


def test_run_locus_coloc_ranks_shared_above_distinct():
    n = 40
    gwas = _gwas_dict(n, peak_idx=20, peak_z=7.0)
    eqtl = {
        "ENSG_SHARED": _eqtl_recs(n, peak_idx=20, peak_z=7.0),   # same peak -> coloc
        "ENSG_DISTINCT": _eqtl_recs(n, peak_idx=3, peak_z=7.0),  # different peak
    }
    res = run_locus_coloc(gwas=gwas, eqtl=eqtl, region="chr1:1-40", min_snps=20)
    assert not res.table.empty
    top = res.table.iloc[0]
    assert top["gene_id"] == "ENSG_SHARED"
    assert top["PP4"] > 0.5
    distinct = res.table[res.table["gene_id"] == "ENSG_DISTINCT"].iloc[0]
    assert distinct["PP4"] < top["PP4"]


def test_run_locus_coloc_min_snps_filter():
    gwas = _gwas_dict(10, peak_idx=5, peak_z=7.0)
    eqtl = {"ENSG_X": _eqtl_recs(10, peak_idx=5, peak_z=7.0)}
    res = run_locus_coloc(gwas=gwas, eqtl=eqtl, region="chr1:1-10", min_snps=20)
    assert res.table.empty  # only 10 shared variants < min_snps=20


def test_run_locus_coloc_allele_flip_matches():
    # eQTL stored with swapped alleles relative to GWAS -> must still match (beta flipped)
    n = 40
    gwas = _gwas_dict(n, peak_idx=20, peak_z=7.0)
    z = np.linspace(-0.5, 0.5, n)
    z[20] = 7.0
    eqtl = {"ENSG_FLIP": [((i + 1, "G", "A"), -z[i] * 0.03, 0.03, 1.0) for i in range(n)]}
    res = run_locus_coloc(gwas=gwas, eqtl=eqtl, region="chr1:1-40", min_snps=20)
    assert not res.table.empty
    assert res.table.iloc[0]["PP4"] > 0.5  # flip handled -> still colocalizes


# ---------------------------------------------------------------------------
# URL helper
# ---------------------------------------------------------------------------


def test_eqtl_catalogue_url_from_qtd():
    url = eqtl_catalogue_url("QTD000251")
    assert url.endswith("QTS000015/QTD000251/QTD000251.all.tsv.gz")


def test_eqtl_catalogue_url_passthrough():
    full = "https://example.org/x.all.tsv.gz"
    assert eqtl_catalogue_url(full) == full


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_gwas_coloc_cli_requires_n_for_finemap(tmp_path):
    from hvantk.tools.qtl.qtlcascade_cli import qtlcascade_group

    r = CliRunner().invoke(
        qtlcascade_group,
        ["gwas-coloc", "--endpoint", "I9_AF", "--chrom", "10", "--lead", "73600000",
         "--eqtl", "QTD000251", "-o", str(tmp_path / "out")],  # fine-map on, no -n
    )
    assert r.exit_code == 1
    # Validation errors go to stderr (err=True); Click 8.1 mixes it into .output,
    # 8.2+ separates it — combine both so the assertion is version-robust.
    combined = r.output
    try:
        combined += r.stderr
    except (ValueError, AttributeError):
        pass
    assert "fine-mapping requires gwas_N and eqtl_N" in combined


def test_gwas_coloc_cli_runs_with_mocked_pipeline(tmp_path, monkeypatch):
    from hvantk.tools.qtl import qtlcascade_cli
    import hvantk.algorithms.qtlcascade.gwas_pipeline as gp

    fake = {
        "region": "chr10:73100000-74100000",
        "gwas": {"min_p_in_region": 2.7e-14},
        "results": {"n_genes_tested": 51, "top_effector": "ENSG00000177791",
                    "top_PP4": 0.812, "report_json": str(tmp_path / "r.json")},
        "fine_map": {"available": True, "credible_sets_gwas": 1,
                     "credible_sets_eqtl": 1, "coloc_susie_PP4": 0.71},
        "verdict": "CONFIRMED (ABF + fine-mapping agree)",
    }
    monkeypatch.setattr(gp, "run_gwas_coloc_pipeline", lambda config: fake)

    r = CliRunner().invoke(
        qtlcascade_cli.qtlcascade_group,
        ["gwas-coloc", "--endpoint", "I9_AF", "--chrom", "10", "--lead", "73600000",
         "--eqtl", "QTD000251", "--no-fine-map", "-o", str(tmp_path)],
    )
    assert r.exit_code == 0, r.output
    assert "CONFIRMED" in r.output
    assert "ENSG00000177791" in r.output


def test_pipeline_user_gene_absent_not_misattributed(tmp_path, monkeypatch):
    # Regression: a user-specified gene absent from the coloc table must NOT
    # inherit the top gene's PP4 (verdict must be about the requested gene).
    import hvantk.algorithms.qtlcascade.gwas_coloc as gc
    from hvantk.algorithms.qtlcascade.gwas_pipeline import (
        GwasColocConfig, run_gwas_coloc_pipeline,
    )

    monkeypatch.setattr(gc, "fetch_finngen_region",
                        lambda *a, **k: _gwas_dict(40, 20, 7.0))
    monkeypatch.setattr(gc, "fetch_eqtl_region",
                        lambda *a, **k: {"ENSG_OTHER": _eqtl_recs(40, 20, 7.0)})
    cfg = GwasColocConfig(
        endpoint="X", chrom="1", lead=1_000_000, eqtl_dataset="QTD",
        gene_of_interest="ENSG_ABSENT", fine_map=False, output_dir=str(tmp_path),
    )
    rep = run_gwas_coloc_pipeline(cfg)
    assert rep["results"]["goi_PP4_abf"] is None          # not the top gene's PP4
    assert rep["results"]["top_effector"] == "ENSG_OTHER"  # top still reported
    assert rep["verdict"].startswith("NO COLOC")
