"""Tests for ``apply_clinvar_training_labels`` / ``ClinVarVariantTableStreamer``."""

from __future__ import annotations

import pytest


@pytest.mark.hail
def test_label_training_set_disease_terms_matches_array_typed_clndn(hail_session):
    """``disease_terms`` must work when ``info.CLNDN`` is ``array<str>``.

    ``hvantk reprocess clinvar:variants`` stores ``info.CLNDN`` as ``array<str>``
    (see tests/snapshots/schema.json), but ``apply_clinvar_training_labels`` used
    to call ``hl.or_else(ht.info.CLNDN, "")`` unconditionally, which raises a Hail
    TypeError (str vs array<str>) the moment ``disease_terms`` is given against a
    table built that way.
    """
    import hail as hl
    from hvantk.skills.clinvar.streamers import apply_clinvar_training_labels

    schema = hl.tstruct(
        id=hl.tstr,
        info=hl.tstruct(
            GENEINFO=hl.tstr,
            CLNSIG=hl.tarray(hl.tstr),
            CLNDN=hl.tarray(hl.tstr),
        ),
    )
    rows = [
        hl.Struct(
            id="matching",
            info=hl.Struct(
                GENEINFO="SOMEGENE:1",
                CLNSIG=["Uncertain_significance"],
                CLNDN=["Hereditary cancer-predisposing syndrome", "not_provided"],
            ),
        ),
        hl.Struct(
            id="non_matching",
            info=hl.Struct(
                GENEINFO="OTHERGENE:2",
                CLNSIG=["Uncertain_significance"],
                CLNDN=["Some_other_condition"],
            ),
        ),
    ]
    ht = hl.Table.parallelize(rows, schema, key="id")

    labeled = apply_clinvar_training_labels(
        ht,
        pathogenic_labels=["Pathogenic"],
        benign_labels=["Benign"],
        disease_terms=["hereditary_cancer-predisposing_syndrome"],
    )
    by_id = {row.id: row for row in labeled.collect()}

    assert by_id["matching"].rf_label == "TP"
    assert "non_matching" not in by_id
