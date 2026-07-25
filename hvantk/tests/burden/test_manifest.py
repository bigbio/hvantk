import yaml

from hvantk.algorithms.burden.manifest import (
    render_cohort_manifest,
    ARCHITECTURE_AXIS_COLUMNS,
)
from hvantk.algorithms.cohort.spec import load_cohort


def test_render_round_trips_via_load_cohort(tmp_path):
    p = tmp_path / "cohort.yaml"
    p.write_text(render_cohort_manifest(name="demo", key="symbol", table="burden.tsv"))
    m = load_cohort(str(p))
    assert m.name == "demo"
    assert m.key == "symbol"
    assert m.key_column == "gene"
    assert m.table == "burden.tsv"
    assert (m.prior.column, m.prior.direction) == ("minp", "lower_is_better")
    assert m.axis("architecture").columns == ARCHITECTURE_AXIS_COLUMNS
    # the declared columns are exactly what `hvantk cohort burden` emits
    assert m.declared_columns() == ("minp",) + ARCHITECTURE_AXIS_COLUMNS


def test_render_honours_key_space_and_key_column(tmp_path):
    txt = render_cohort_manifest(
        name="d", key="gene_id", table="t.tsv", key_column="ensembl"
    )
    doc = yaml.safe_load(txt)  # the leading `#` header lines are YAML comments
    assert doc["key"] == "gene_id"
    assert doc["key_column"] == "ensembl"
    assert doc["prior"] == {"column": "minp", "direction": "lower_is_better"}


def test_prior_column_is_overridable():
    """A consumed prior is rarely called 'minp'.

    `hvantk cohort burden` names its own output column `minp`, but the manifest contract
    also describes cohorts whose prior table was produced elsewhere -- Epi25 DEE calls it
    `p_dee`, SCHEMA calls it `p_burden`. Hardcoding `minp` emits a manifest that
    validates and then resolves to a column that does not exist.
    """
    import yaml

    from hvantk.algorithms.burden.manifest import render_cohort_manifest

    doc = yaml.safe_load(
        render_cohort_manifest(
            name="epilepsy_dee", key="symbol", table="/x/epi25.tsv",
            prior_column="p_dee",
        )
    )
    assert doc["prior"]["column"] == "p_dee"
    assert doc["prior"]["direction"] == "lower_is_better"


def test_architecture_axis_can_be_omitted():
    """A gene-level-only cohort has no architecture columns to point at.

    Epi25 DEE ships per-gene p-values and nothing else -- no cohort variant table, so
    no n_case_var/conc/driver_af. Emitting the architecture axis anyway produces a
    manifest that names columns the table does not have, and the absence of an artifact
    veto for that cohort is a real property that must be visible, not papered over.
    """
    import yaml

    from hvantk.algorithms.burden.manifest import render_cohort_manifest

    doc = yaml.safe_load(
        render_cohort_manifest(
            name="epilepsy_dee", key="symbol", table="/x/epi25.tsv",
            prior_column="p_dee", architecture_columns=None,
        )
    )
    assert "cohort_axes" not in doc


def test_defaults_are_unchanged():
    """The burden op's own emission must be byte-identical to before."""
    import yaml

    from hvantk.algorithms.burden.manifest import (
        ARCHITECTURE_AXIS_COLUMNS,
        render_cohort_manifest,
    )

    doc = yaml.safe_load(
        render_cohort_manifest(name="chd", key="symbol", table="/x/chd.tsv")
    )
    assert doc["prior"]["column"] == "minp"
    assert doc["cohort_axes"] == [
        {"axis": "architecture", "columns": list(ARCHITECTURE_AXIS_COLUMNS)}
    ]
