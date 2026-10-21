"""Tests for ``apply_clinvar_training_labels`` / ``ClinVarVariantTableStreamer``."""

from __future__ import annotations

import pytest


@pytest.mark.hail
def test_label_training_set_disease_terms_matches_array_typed_clndn(hail_session):
    """``disease_terms`` must match whole CLNDN names, not comma fragments.

    ``hvantk reprocess clinvar:variants`` stores ``info.CLNDN`` as ``array<str>``
    (see tests/snapshots/schema.json) because the VCF header declares CLNDN
    ``Number=.``: ``hl.import_vcf`` splits the raw value at every comma, including
    a comma inside a single disease name. ``CLNDN=Glycogen_storage_disease,_type_IV
    |Polyglucosan_body_myopathy_type_1`` therefore arrives as
    ``["Glycogen_storage_disease", "_type_IV|Polyglucosan_body_myopathy_type_1"]`` --
    comma fragments of one ``'|'``-separated list, not independent whole names. The
    full name must still match, and the bare fragment ("Glycogen storage disease")
    must not.
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
            id="comma_variant",
            info=hl.Struct(
                GENEINFO="SOMEGENE:1",
                CLNSIG=["Uncertain_significance"],
                CLNDN=[
                    "Glycogen_storage_disease",
                    "_type_IV|Polyglucosan_body_myopathy_type_1",
                ],
            ),
        ),
    ]
    ht = hl.Table.parallelize(rows, schema, key="id")

    # The whole, comma-containing name matches.
    full_name = apply_clinvar_training_labels(
        ht,
        pathogenic_labels=["Pathogenic"],
        benign_labels=["Benign"],
        disease_terms=["Glycogen storage disease, type IV"],
    )
    assert {row.id for row in full_name.collect()} == {"comma_variant"}

    # The comma fragment alone (a false positive under the old splitting) must not.
    fragment_only = apply_clinvar_training_labels(
        ht,
        pathogenic_labels=["Pathogenic"],
        benign_labels=["Benign"],
        disease_terms=["Glycogen storage disease"],
    )
    assert fragment_only.count() == 0
