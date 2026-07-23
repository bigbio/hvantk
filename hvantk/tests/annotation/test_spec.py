"""Feature-spec parsing and validation (P2a: gene_id-keyed entries only)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

SCHEMA = Path("hvantk/resources/schemas/feature_spec.schema.json")


# These tests validate through jsonschema's version-tolerant entry points
# (``validator_for`` + ``validate``) rather than the 4.x-only ``Draft202012Validator``
# class. CI's conda environment pins an older jsonschema (3.2.0) that has no
# ``Draft202012Validator``, and the production loader (``spec.load_spec``) and every existing
# plugin/tool loader already validate this way, so the tests must match that path. The
# assertions here (required fields, enum, additionalProperties) are all draft-4+ features, so
# they hold whether the schema is validated as 2020-12 or under an older draft's fallback.


def test_feature_spec_schema_is_valid():
    import jsonschema

    schema = json.loads(SCHEMA.read_text())
    assert schema["$schema"] == "https://json-schema.org/draft/2020-12/schema"
    # Must not raise: the schema is a well-formed schema for whatever draft this jsonschema
    # resolves it to (validator_for picks 2020-12 on 4.x, an older fallback on 3.x).
    validator_cls = jsonschema.validators.validator_for(schema)
    validator_cls.check_schema(schema)


def test_schema_accepts_a_gene_id_entry_and_rejects_a_bad_one():
    import jsonschema
    import pytest

    schema = json.loads(SCHEMA.read_text())

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
    jsonschema.validate(good, schema)  # must not raise

    bad = {
        "name": "x",
        "layer1": [{"axis": "constraint"}],
    }  # missing source/key/columns
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


def test_schema_rejects_a_non_gene_id_key():
    """P2a restricts key to gene_id; a well-formed hgnc_id entry must fail validation.

    This is the schema half of the key-restriction defense-in-depth. P2c widens the enum;
    when it does, this test is the deliberate tripwire it must update.
    """
    import jsonschema
    import pytest

    schema = json.loads(SCHEMA.read_text())
    entry = {
        "name": "x",
        "layer1": [
            {
                "axis": "gene-disease",
                "source": "clingen:gene-disease",
                "key": "hgnc_id",
                "columns": ["classification"],
            }
        ],
    }
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(entry, schema)


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

    with pytest.raises(KeyError):
        load_spec(p).entry("expression")


def test_a_multi_entry_spec_addresses_each_axis_independently(tmp_path):
    from hvantk.algorithms.annotation.spec import load_spec

    doc = (
        "name: t\n"
        "layer1:\n"
        "  - {axis: constraint, source: gnomad-metrics:metrics, key: gene_id, "
        "columns: [mis_z]}\n"
        "  - {axis: gevir, source: gevir:metrics, key: gene_id, "
        "columns: [gevir_pct, virlof_pct]}\n"
    )
    p = tmp_path / "s.yaml"
    p.write_text(doc)

    spec = load_spec(p)
    assert spec.entry("constraint").source == "gnomad-metrics:metrics"
    gevir = spec.entry("gevir")
    assert gevir.source == "gevir:metrics"
    assert gevir.columns == ("gevir_pct", "virlof_pct")


def test_duplicate_axis_labels_are_rejected(tmp_path):
    from hvantk.algorithms.annotation.spec import load_spec

    doc = (
        "name: t\n"
        "layer1:\n"
        "  - {axis: constraint, source: a:x, key: gene_id, columns: [c1]}\n"
        "  - {axis: constraint, source: b:y, key: gene_id, columns: [c2]}\n"
    )
    p = tmp_path / "s.yaml"
    p.write_text(doc)

    with pytest.raises(ValueError, match="duplicate axis"):
        load_spec(p)
