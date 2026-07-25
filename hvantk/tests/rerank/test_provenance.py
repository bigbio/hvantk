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


def test_unrecognised_source_conflicts_with_matching_unrecognised_source():
    """The other half of self-mapping: an unrecognised name must still conflict with itself.

    ``test_different_classes_do_not_conflict`` only proves an unmapped source does not
    conflict with an UNRELATED class. That alone is also satisfied by a resolver that
    silently drops any source with no declared equivalence class instead of mapping it
    onto itself -- dropping never conflicts with anything, from either side. Pin the other
    half: a feature trained on an unmapped source and a label derived from that SAME
    unmapped source must still be reported as conflicted.
    """
    a = resolve_arms(
        {"toy_score": frozenset({"simulated"})}, frozenset({"simulated"}), DEFAULT_EQUIVALENCE
    )
    assert a.clean == ()
    assert a.conflicted == {"toy_score": frozenset({"simulated"})}


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


def test_custom_equivalence_map_changes_the_verdict():
    """``equivalence`` is live configuration, not a comment about DEFAULT_EQUIVALENCE.

    Every other test in this file passes DEFAULT_EQUIVALENCE, so an implementation that
    ignored the parameter and closed over the module constant instead would still pass all
    of them. Prove the parameter is actually load-bearing: a custom map groups two sources
    that DEFAULT_EQUIVALENCE treats as unrelated, and the SAME raw inputs must flip from
    clean to conflicted purely because of which equivalence map was supplied -- this is
    what makes "extending the vocabulary is a data edit, not a code change" a real property.
    """
    custom_equivalence = {"lab_x_predictor": ["FooScore", "BarLabel"]}
    conflicted = resolve_arms(
        {"widget": frozenset({"FooScore"})}, frozenset({"BarLabel"}), custom_equivalence
    )
    assert conflicted.clean == ()
    assert conflicted.conflicted == {"widget": frozenset({"lab_x_predictor"})}

    # Identical raw inputs, only the equivalence map differs: DEFAULT_EQUIVALENCE has no
    # class containing "FooScore" or "BarLabel", so each maps to itself and they don't
    # collide -- the verdict flips to clean.
    clean = resolve_arms(
        {"widget": frozenset({"FooScore"})}, frozenset({"BarLabel"}), DEFAULT_EQUIVALENCE
    )
    assert clean.clean == ("widget",)
