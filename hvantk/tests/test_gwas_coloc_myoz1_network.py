"""End-to-end positive-control regression: AF 10q22 → MYOZ1 must be CONFIRMED.

Network only (FinnGen + eQTL Catalogue + 1000G LD, all via remote tabix). The
fine-mapping is now pure-Python SuSiE-RSS + coloc.susie (no R/bcftools). Auto-marked
``network`` by the filename suffix, so it is excluded from the default fast
run; execute with ``pytest -m network hvantk/tests/test_gwas_coloc_myoz1_network.py``.

This guards the whole chain: ABF coloc recovers MYOZ1 as the top effector AND
SuSiE/coloc.susie confirms it — the discrimination that makes the workflow
trustworthy (cf. the CHD 17q21/NSF look-alike, which is REFUTED; see
``test_gwas_coloc_nsf_network.py``).
"""

import pytest

from hvantk.algorithms.qtlcascade.finemap import finemap_available
from hvantk.algorithms.qtlcascade.gwas_pipeline import (
    GwasColocConfig,
    run_gwas_coloc_pipeline,
)

MYOZ1 = "ENSG00000177791"


@pytest.mark.skipif(not finemap_available()[0],
                    reason="needs pysam for the 1000G LD reference")
def test_af_myoz1_positive_control_confirmed(tmp_path):
    cfg = GwasColocConfig(
        endpoint="I9_AF", chrom="10", lead=73600000,
        eqtl_dataset="QTD000251",            # GTEx heart atrial appendage
        gene_of_interest=MYOZ1,
        gwas_N=261395, eqtl_N=372, fine_map=True,
        output_dir=str(tmp_path), ld_cache_dir=str(tmp_path / "ld"),
    )
    report = run_gwas_coloc_pipeline(cfg)

    # ABF: MYOZ1 is the top effector with high H4
    assert report["results"]["top_effector"] == MYOZ1
    assert report["results"]["top_PP4"] > 0.5

    # Fine-mapping confirms it (credible set in each trait + shared)
    fm = report["fine_map"]
    assert fm["available"] is True
    assert fm["credible_sets_gwas"] >= 1
    assert fm["credible_sets_eqtl"] >= 1
    assert fm["coloc_susie_PP4"] is not None and fm["coloc_susie_PP4"] > 0.5

    assert report["verdict"].startswith("CONFIRMED")
