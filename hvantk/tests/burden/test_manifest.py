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
