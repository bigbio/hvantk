import pytest

from hvantk.algorithms.rerank.provenance import DEFAULT_EQUIVALENCE, resolve_arms


def test_empty_trained_on_is_clean():
    """A predictor trained on nothing label-derived cannot be circular."""
    a = resolve_arms({"phyloP": frozenset()}, frozenset({"GenCC"}), DEFAULT_EQUIVALENCE)
    assert a.clean == ("phyloP",)
    assert a.conflicted == {}
    assert a.undeclared == ()


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


def test_undeclared_is_conflicted_and_not_clean():
    """Absent provenance must fail safe: usable, but never in the headline arm.

    The bar is recorded as a conflict rather than as a third neutral-sounding bucket. A
    column nobody has vouched for is presumed to have seen whatever the label derives
    from, so the presumed conflict is the label's own classes.
    """
    a = resolve_arms({"mystery": None}, frozenset({"GenCC"}), DEFAULT_EQUIVALENCE)
    assert a.clean == ()
    assert a.conflicted == {"mystery": frozenset({"curated_disease_db"})}
    assert a.undeclared == ("mystery",)


def test_undeclared_columns_are_nameable_for_review():
    """`undeclared` is a subset of `conflicted`, not a parallel bucket.

    Without the list an undeclared column is indistinguishable from a genuinely circular
    one, and nobody can tell which manifest entry is missing -- which is how four dbNSFP
    columns (MutFormer, PHACTboost) sat unnoticed outside the clean arm.
    """
    a = resolve_arms(
        {"REVEL": frozenset({"ClinVar"}), "mystery": None},
        frozenset({"GenCC"}),
        DEFAULT_EQUIVALENCE,
    )
    assert set(a.undeclared) <= set(a.conflicted)
    assert "REVEL" in a.conflicted and "REVEL" not in a.undeclared


def test_all_columns_does_not_double_count_undeclared():
    """`undeclared` members are already `conflicted` keys; counting both duplicates them."""
    a = resolve_arms(
        {"phyloP": frozenset(), "mystery": None},
        frozenset({"GenCC"}),
        DEFAULT_EQUIVALENCE,
    )
    assert sorted(a.all_columns) == ["mystery", "phyloP"]
    assert len(a.all_columns) == len(set(a.all_columns))


def test_all_columns_is_clean_plus_conflicted():
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


def test_dbnsfp_declares_provenance_for_its_predictors():
    """The plugin author declares training sources once; every consumer inherits them."""
    from pathlib import Path

    import yaml

    import hvantk

    manifest = Path(hvantk.__file__).parent / "skills" / "dbnsfp" / "plugin.yaml"
    doc = yaml.safe_load(manifest.read_text())
    scores = next(d for d in doc["datasets"] if d["name"] == "variants")["scores"]

    # disease-database trained -> must be declared, else they silently reach the headline
    assert set(scores["REVEL_rankscore"]["trained_on"]) >= {"HGMD"}
    assert "ClinVar" in scores["MVP_rankscore"]["trained_on"]
    # trained on simulated alleles, not curated disease sets
    assert scores["CADD_raw_rankscore"]["trained_on"] == ["simulated"]
    # pure conservation: nothing to conflict with
    assert scores["phyloP100way_vertebrate_rankscore"]["trained_on"] == []


def test_declared_dbnsfp_scores_resolve_into_arms_against_a_clinvar_label():
    """End-to-end: the declarations actually drive the clean/all split.

    Pins the join between the plugin manifest and `resolve_arms` -- a declaration that
    parses but never reaches the resolver would be documentation, not a control.
    """
    from pathlib import Path

    import yaml

    import hvantk
    from hvantk.algorithms.rerank.provenance import DEFAULT_EQUIVALENCE, resolve_arms

    manifest = Path(hvantk.__file__).parent / "skills" / "dbnsfp" / "plugin.yaml"
    doc = yaml.safe_load(manifest.read_text())
    scores = next(d for d in doc["datasets"] if d["name"] == "variants")["scores"]
    provenance = {k: frozenset(v["trained_on"]) for k, v in scores.items()}

    arms = resolve_arms(provenance, frozenset({"ClinGen"}), DEFAULT_EQUIVALENCE)

    assert "REVEL_rankscore" in arms.conflicted          # HGMD/ClinVar vs a ClinGen label
    assert "phyloP100way_vertebrate_rankscore" in arms.clean
    assert "CADD_raw_rankscore" in arms.clean            # 'simulated' is not a disease db
    assert arms.undeclared == ()                         # every declared score is declared


def test_a_source_in_two_equivalence_classes_is_rejected():
    """Ambiguous vocabulary must fail loud, not resolve by dict iteration order.

    A duplicate member would otherwise be assigned to whichever class the dict happened
    to process last, silently moving every feature trained on that source between the
    clean and all arms. Same failure shape the CLI rejects for a repeated `--prepared`
    axis (PR #229) -- reject the ambiguity rather than guess.
    """
    bad = {"curated_disease_db": ["ClinVar", "HGMD"], "other": ["ClinVar"]}
    with pytest.raises(ValueError, match="two classes"):
        resolve_arms({"x": frozenset({"ClinVar"})}, frozenset({"HGMD"}), bad)


def test_a_source_repeated_within_one_class_is_fine():
    """Only a CROSS-class duplicate is ambiguous; repeating inside one class is not."""
    dup = {"curated_disease_db": ["ClinVar", "ClinVar", "HGMD"]}
    arms = resolve_arms({"x": frozenset({"ClinVar"})}, frozenset({"HGMD"}), dup)
    assert "x" in arms.conflicted
