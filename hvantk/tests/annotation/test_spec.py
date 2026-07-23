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
