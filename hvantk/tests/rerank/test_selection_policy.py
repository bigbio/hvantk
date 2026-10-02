import pytest


def test_empty_document_is_a_valid_policy(tmp_path):
    """An empty selection.yaml must yield working defaults, not an error."""
    from hvantk.algorithms.rerank.provenance import DEFAULT_EQUIVALENCE
    from hvantk.algorithms.rerank.selection import load_policy

    p = tmp_path / "selection.yaml"
    p.write_text("{}\n")
    policy, equivalence = load_policy(p)
    assert policy.q == 0.10 and policy.wrapper == "none" and policy.inner_folds == 3
    assert equivalence == DEFAULT_EQUIVALENCE


def test_overrides_are_applied(tmp_path):
    from hvantk.algorithms.rerank.selection import load_policy

    p = tmp_path / "selection.yaml"
    p.write_text("q: 0.05\nwrapper: rfecv\ninner_folds: 5\n")
    policy, _ = load_policy(p)
    assert policy.q == 0.05 and policy.wrapper == "rfecv" and policy.inner_folds == 5


def test_equivalence_map_is_overridable(tmp_path):
    """Extending the vocabulary must be a data edit, not a code change."""
    from hvantk.algorithms.rerank.selection import load_policy

    p = tmp_path / "selection.yaml"
    p.write_text(
        "equivalence:\n  curated_disease_db: [ClinVar, GenCC]\n  my_class: [SourceA]\n"
    )
    _, equivalence = load_policy(p)
    assert equivalence["my_class"] == ["SourceA"]


def test_unknown_key_is_rejected(tmp_path):
    import jsonschema

    from hvantk.algorithms.rerank.selection import load_policy

    p = tmp_path / "selection.yaml"
    p.write_text("qq: 0.05\n")
    with pytest.raises(jsonschema.ValidationError):
        load_policy(p)


def test_out_of_range_q_is_rejected(tmp_path):
    import jsonschema

    from hvantk.algorithms.rerank.selection import load_policy

    p = tmp_path / "selection.yaml"
    p.write_text("q: 1.5\n")
    with pytest.raises(jsonschema.ValidationError):
        load_policy(p)


def test_equivalence_is_not_passed_to_the_policy_dataclass(tmp_path):
    """`equivalence` is provenance vocabulary, not a filter knob.

    Pinned because the loader pops it out of the document before constructing the
    policy; if that pop is ever lost the dataclass raises an opaque TypeError instead
    of a schema error, which is a confusing failure for a valid file.
    """
    from hvantk.algorithms.rerank.selection import load_policy

    p = tmp_path / "selection.yaml"
    p.write_text("q: 0.05\nequivalence:\n  my_class: [SourceA]\n")
    policy, equivalence = load_policy(p)
    assert policy.q == 0.05
    assert not hasattr(policy, "equivalence")
    assert equivalence == {"my_class": ["SourceA"]}


# --- the method vocabulary -------------------------------------------------------------

_METHOD_FIELDS = ("univariate", "redundancy", "wrapper", "wrapper_estimator")


def _schema_enums():
    import json

    from hvantk.algorithms.rerank.selection import _POLICY_SCHEMA_PATH

    props = json.loads(_POLICY_SCHEMA_PATH.read_text())["properties"]
    return {f: tuple(props[f]["enum"]) for f in _METHOD_FIELDS}


def test_method_vocabularies_are_the_schema_enums():
    """ONE vocabulary: the module constants, pinned to the schema's enums.

    `load_policy` validates YAML against the schema and the dataclass validates the
    Python constructor against the constants; if the two ever disagree, one path accepts
    a spelling the other rejects, and a config that works from YAML stops working from
    Python (or the reverse) with nothing to say why.
    """
    from hvantk.algorithms.rerank import selection as s

    assert _schema_enums() == {
        "univariate": s.UNIVARIATE_METHODS,
        "redundancy": s.REDUNDANCY_METHODS,
        "wrapper": s.WRAPPER_METHODS,
        "wrapper_estimator": s.WRAPPER_ESTIMATORS,
    }


@pytest.mark.parametrize(
    "field,value",
    [
        ("univariate", "AUC"),
        ("univariate", "pearson"),
        ("redundancy", "pearson"),
        ("redundancy", "Spearman"),
        ("wrapper", "RFECV"),
        ("wrapper", "rfe"),
        ("wrapper", None),
        ("wrapper_estimator", "svm"),
    ],
)
def test_a_misspelled_method_is_rejected_not_silently_off(field, value):
    """`select_axis` branches on exact string equality with no else, so before this check
    `wrapper="RFECV"` or `redundancy="pearson"` silently took the "off" branch: the step
    the caller asked for never ran, every column survived it, and nothing said so (two
    noise columns were "selected" under three typos where the correct spelling dropped
    both). The error must name the field, the value and the allowed spellings."""
    from hvantk.algorithms.rerank.selection import SelectionPolicy

    with pytest.raises(ValueError) as info:
        SelectionPolicy(**{field: value})
    msg = str(info.value)
    assert f"SelectionPolicy.{field}" in msg and repr(value) in msg
    for allowed in _schema_enums()[field]:
        assert repr(allowed) in msg


def test_every_allowed_spelling_constructs():
    from hvantk.algorithms.rerank.selection import SelectionPolicy

    SelectionPolicy()  # the defaults are themselves allowed spellings
    for field, allowed in _schema_enums().items():
        for value in allowed:
            assert getattr(SelectionPolicy(**{field: value}), field) == value


def test_a_misspelled_method_in_yaml_is_a_schema_error_and_the_loader_round_trips(
    tmp_path,
):
    """The YAML path rejects the same typo first, as a schema error, and still round-trips
    every spelling the constructor accepts."""
    import jsonschema

    from hvantk.algorithms.rerank.selection import load_policy

    p = tmp_path / "selection.yaml"
    p.write_text("wrapper: RFECV\n")
    with pytest.raises(jsonschema.ValidationError):
        load_policy(p)

    p.write_text(
        "univariate: none\nredundancy: none\nwrapper: rfecv\n"
        "wrapper_estimator: random_forest\n"
    )
    policy, _ = load_policy(p)
    assert (policy.univariate, policy.redundancy, policy.wrapper) == (
        "none",
        "none",
        "rfecv",
    )
    assert policy.wrapper_estimator == "random_forest"


# --- config wiring ---------------------------------------------------------------------


def test_config_defaults_to_no_selection():
    from hvantk.algorithms.rerank.config import Config

    assert Config.__dataclass_fields__["selection"].default is None


@pytest.mark.parametrize("bad", [42, "auc", {"wrapper": "rfecv"}])
def test_config_rejects_a_non_selection_policy(bad):
    """A bare int, string or dict here used to construct and then fail deep inside a
    per-fold selector, after scoring had begun, with an AttributeError that named none of
    this -- the failure mode Config.leakage's, Config.nulls's and Config.blocks's type
    checks exist for."""
    from hvantk.algorithms.rerank.config import Config

    with pytest.raises(TypeError, match="SelectionPolicy"):
        Config(name="x", features=[], labels=None, selection=bad).__post_init__()
