import pytest

from hvantk.algorithms.rerank.provenance import DEFAULT_EQUIVALENCE, resolve_arms


def test_empty_trained_on_is_clean():
    """A predictor trained on nothing label-derived cannot be circular."""
    a = resolve_arms({"phyloP": frozenset()}, frozenset({"GenCC"}), DEFAULT_EQUIVALENCE)
    assert a.clean == ("phyloP",)
    assert a.conflicted == {}
    assert a.unknown == ()


def test_shared_equivalence_class_conflicts():
    """ClinVar-trained vs GenCC-derived: different names, same curation class."""
    a = resolve_arms(
        {"REVEL": frozenset({"HGMD", "ClinVar"})}, frozenset({"GenCC"}), DEFAULT_EQUIVALENCE
    )
    assert a.clean == ()
    assert "REVEL" in a.conflicted
    assert a.conflicted["REVEL"] == frozenset({"curated_disease_db"})


def test_different_classes_do_not_conflict():
    """CADD is trained on simulated alleles; no curation class is involved."""
    a = resolve_arms(
        {"CADD": frozenset({"simulated"})}, frozenset({"GenCC"}), DEFAULT_EQUIVALENCE
    )
    assert a.clean == ("CADD",)


def test_undeclared_is_unknown_and_not_clean():
    """Absent provenance must fail safe: usable, but never in the headline arm."""
    a = resolve_arms({"mystery": None}, frozenset({"GenCC"}), DEFAULT_EQUIVALENCE)
    assert a.clean == ()
    assert a.unknown == ("mystery",)
    assert a.conflicted == {}


def test_all_columns_is_clean_plus_conflicted_plus_unknown():
    a = resolve_arms(
        {"phyloP": frozenset(), "REVEL": frozenset({"ClinVar"}), "mystery": None},
        frozenset({"GenCC"}),
        DEFAULT_EQUIVALENCE,
    )
    assert set(a.all_columns) == {"phyloP", "REVEL", "mystery"}
    assert a.clean == ("phyloP",)


def test_label_with_no_declared_provenance_conflicts_with_nothing():
    """A burden-only label shares no class with a disease-trained predictor."""
    a = resolve_arms(
        {"REVEL": frozenset({"ClinVar"})}, frozenset(), DEFAULT_EQUIVALENCE
    )
    assert a.clean == ("REVEL",)
