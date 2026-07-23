"""Feature-spec parsing and validation (P2a: gene_id-keyed entries only)."""
from __future__ import annotations

import json
from pathlib import Path

SCHEMA = Path("hvantk/resources/schemas/feature_spec.schema.json")


def test_feature_spec_schema_is_valid_draft_2020_12():
    import jsonschema

    schema = json.loads(SCHEMA.read_text())
    assert schema["$schema"] == "https://json-schema.org/draft/2020-12/schema"
    # Must not raise: the schema itself is a well-formed schema.
    jsonschema.Draft202012Validator.check_schema(schema)


def test_schema_accepts_a_gene_id_entry_and_rejects_a_bad_one():
    import jsonschema

    schema = json.loads(SCHEMA.read_text())
    validator = jsonschema.Draft202012Validator(schema)

    good = {
        "name": "chd-v2",
        "layer1": [
            {
                "axis": "constraint",
                "source": "gnomad-metrics:metrics",
                "key": "gene_id",
                "columns": ["mis_z", "pLI"],
                "min_mapping_rate": 0.9,
            }
        ],
    }
    assert validator.is_valid(good), list(validator.iter_errors(good))

    bad = {
        "name": "x",
        "layer1": [{"axis": "constraint"}],
    }  # missing source/key/columns
    assert not validator.is_valid(bad)


def test_load_spec_parses_a_gene_id_entry(tmp_path):
    from hvantk.algorithms.annotation.spec import load_spec

    p = tmp_path / "s.yaml"
    p.write_text(
        "name: chd-v2\n"
        "layer1:\n"
        "  - axis: constraint\n"
        "    source: gnomad-metrics:metrics\n"
        "    key: gene_id\n"
        "    columns: [mis_z, pLI]\n"
        "    min_mapping_rate: 0.9\n"
    )
    spec = load_spec(p)
    assert spec.name == "chd-v2"
    entry = spec.entry("constraint")
    assert entry.source == "gnomad-metrics:metrics"
    assert entry.key == "gene_id"
    assert entry.columns == ("mis_z", "pLI")
    assert entry.min_mapping_rate == 0.9


def test_min_mapping_rate_defaults_to_0_9(tmp_path):
    from hvantk.algorithms.annotation.spec import load_spec

    p = tmp_path / "s.yaml"
    p.write_text(
        "name: x\n"
        "layer1:\n"
        "  - axis: constraint\n"
        "    source: gnomad-metrics:metrics\n"
        "    key: gene_id\n"
        "    columns: [mis_z]\n"
    )
    assert load_spec(p).entry("constraint").min_mapping_rate == 0.9


def test_load_spec_rejects_an_invalid_spec(tmp_path):
    import jsonschema

    from hvantk.algorithms.annotation.spec import load_spec

    p = tmp_path / "bad.yaml"
    p.write_text(
        "name: x\nlayer1:\n  - axis: constraint\n"
    )  # missing source/key/columns
    import pytest

    with pytest.raises(jsonschema.ValidationError):
        load_spec(p)


def test_entry_raises_for_a_missing_axis(tmp_path):
    from hvantk.algorithms.annotation.spec import load_spec

    p = tmp_path / "s.yaml"
    p.write_text(
        "name: x\n"
        "layer1:\n"
        "  - axis: constraint\n"
        "    source: gnomad-metrics:metrics\n"
        "    key: gene_id\n"
        "    columns: [mis_z]\n"
    )
    import pytest

    with pytest.raises(KeyError):
        load_spec(p).entry("expression")
