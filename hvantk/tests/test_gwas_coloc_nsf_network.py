"""End-to-end negative-control regression: CHD 17q21 → NSF must be REFUTED.

The septal-defect 17q21 locus has an *identical-looking* single-variant ABF
signal (PP4 ≈ 0.81 for NSF) that does NOT survive SuSiE-RSS + coloc.susie
fine-mapping — neither trait yields a credible set, so coloc.susie PP4 = 0 and
the verdict is REFUTED. This is the discrimination that makes the fine-mapping
layer worth running (cf. the AF→MYOZ1 positive control, which is CONFIRMED).

Network only (FinnGen + eQTL Catalogue + 1000G LD via remote tabix); pure-Python
fine-mapping (no R/bcftools). Auto-marked ``network`` by the filename suffix;
run with ``pytest -m network hvantk/tests/test_gwas_coloc_nsf_network.py``.
"""

import pytest

from hvantk.algorithms.qtlcascade.finemap import finemap_available
from hvantk.algorithms.qtlcascade.gwas_pipeline import (
    GwasColocConfig,
    run_gwas_coloc_pipeline,
)

NSF = "ENSG00000073969"


@pytest.mark.skipif(not finemap_available()[0],
                    reason="needs pysam for the 1000G LD reference")
def test_chd_nsf_look_alike_refuted(tmp_path):
    cfg = GwasColocConfig(
        endpoint="Q17_SEPTA_DEFEC", chrom="17", lead=46890164,
        eqtl_dataset="QTD000136",            # GTEx (eQTL Catalogue) cis-eQTL
        gene_of_interest=NSF,
        gwas_N=412181, eqtl_N=213, fine_map=True,
        output_dir=str(tmp_path), ld_cache_dir=str(tmp_path / "ld"),
    )
    report = run_gwas_coloc_pipeline(cfg)

    # ABF over-calls: NSF looks like the top effector with high H4 …
    assert report["results"]["top_effector"] == NSF
    assert report["results"]["top_PP4"] > 0.5

    # … but fine-mapping finds no credible set in either trait -> REFUTED.
    fm = report["fine_map"]
    assert fm["available"] is True
    assert fm["credible_sets_gwas"] == 0
    assert fm["credible_sets_eqtl"] == 0
    assert fm["coloc_susie_PP4"] == 0.0

    assert report["verdict"].startswith("REFUTED")
