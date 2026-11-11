"""Tests for ``label_clinvar_variants`` / ``ClinVarVariantTableStreamer``."""

from __future__ import annotations

import pytest


@pytest.mark.hail
def test_review_stars_rejoin_split_clnrevstat(hail_session):
    """Review stars must be read from the whole CLNREVSTAT, not its comma fragments.

    The VCF header declares CLNREVSTAT ``Number=.``, so ``hl.import_vcf`` splits it at
    every comma: ``criteria_provided,_single_submitter`` arrives as
    ``["criteria_provided", "_single_submitter"]``. Looked up fragment by fragment, no
    fragment is a known review status and every row would get 0 stars -- dropping all
    rows under a ">= 1 star" filter. The fragments must be re-joined with ``,`` first.
    """
    import hail as hl

    from hvantk.skills.clinvar.streamers import label_clinvar_variants

    schema = hl.tstruct(
        locus=hl.tlocus("GRCh38"),
        alleles=hl.tarray(hl.tstr),
        info=hl.tstruct(
            CLNSIG=hl.tarray(hl.tstr),
            CLNREVSTAT=hl.tarray(hl.tstr),
            MC=hl.tarray(hl.tstr),
            GENEINFO=hl.tstr,
        ),
    )
    rows = [
        hl.Struct(
            locus=hl.Locus("chr1", 100000, "GRCh38"),
            alleles=["C", "T"],
            info=hl.Struct(
                CLNSIG=["Pathogenic"],
                CLNREVSTAT=["criteria_provided", "_single_submitter"],
                MC=["SO:0001583|missense_variant"],
                GENEINFO="SOMEGENE:1",
            ),
        )
    ]
    ht = hl.Table.parallelize(rows, schema, key=["locus", "alleles"])

    (row,) = label_clinvar_variants(
        ht, pathogenic=["Pathogenic"], benign=["Benign"]
    ).collect()

    assert row.clinvar_review_status == "criteria_provided,_single_submitter"
    assert row.clinvar_review_stars == 1
