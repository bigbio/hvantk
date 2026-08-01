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
    p.write_text("equivalence:\n  curated_disease_db: [ClinVar, GenCC]\n  my_class: [SourceA]\n")
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
