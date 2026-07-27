"""Feature-spec parsing and validation (keys: gene_id, hgnc_id, symbol, variant)."""
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


def test_schema_rejects_an_unknown_key():
    """Validates that only the allowed keys (gene_id, hgnc_id, symbol) pass validation.

    This is the schema half of the key-restriction defense-in-depth. P2c widened the enum
    to gene_id/hgnc_id/symbol; this tripwire now asserts a still-unknown key (protein_id) is
    rejected, and must be re-pointed at a new unknown key if the enum widens again.
    """
    import jsonschema
    import pytest

    schema = json.loads(SCHEMA.read_text())
    entry = {
        "name": "x",
        "layer1": [
            {
                "axis": "protein",
                "source": "some:source",
                "key": "protein_id",
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


def test_hgnc_id_and_symbol_keys_are_accepted(tmp_path):
    from hvantk.algorithms.annotation.spec import load_spec

    doc = (
        "name: t\n"
        "layer1:\n"
        "  - {axis: a, source: s:one, key: hgnc_id, columns: [c1]}\n"
        "  - {axis: b, source: s:two, key: symbol, columns: [c2]}\n"
    )
    p = tmp_path / "s.yaml"
    p.write_text(doc)

    spec = load_spec(p)
    assert spec.entry("a").key == "hgnc_id"
    assert spec.entry("b").key == "symbol"


def test_an_unknown_key_is_still_rejected(tmp_path):
    import jsonschema

    from hvantk.algorithms.annotation.spec import load_spec

    doc = (
        "name: t\n"
        "layer1:\n"
        "  - {axis: a, source: s:one, key: protein_id, columns: [c1]}\n"
    )
    p = tmp_path / "s.yaml"
    p.write_text(doc)

    with pytest.raises(jsonschema.ValidationError):
        load_spec(p)


def test_variant_key_with_aggregate_parses(tmp_path):
    from hvantk.algorithms.annotation.spec import load_spec

    doc = (
        "name: t\n"
        "layer1:\n"
        "  - axis: tolerance\n"
        "    source: dbnsfp:variants\n"
        "    key: variant\n"
        "    columns: [revel_mean, revel_frac_gt_0.5]\n"
        "    aggregate:\n"
        "      by: Ensembl_geneid\n"
        "      to: gene_id\n"
        "      filter: missense\n"
        "      scores:\n"
        "        revel: {column: REVEL_score, stats: [mean, frac_gt_0.5]}\n"
    )
    p = tmp_path / "s.yaml"
    p.write_text(doc)
    spec = load_spec(p)
    e = spec.entry("tolerance")
    assert e.key == "variant"
    assert e.aggregate.by == "Ensembl_geneid" and e.aggregate.to == "gene_id"
    assert e.aggregate.reduce == "max"  # default
    assert e.aggregate.scores[0].name == "revel"
    assert e.aggregate.scores[0].column == "REVEL_score"
    assert e.aggregate.scores[0].stats == ("mean", "frac_gt_0.5")
    assert e.aggregate.count_name == "n_possible_missense"  # dbNSFP-era default


def test_aggregate_count_name_is_configurable(tmp_path):
    """Non-dbNSFP sources name their own row-count column (PTM sites, eQTL pairs)."""
    from hvantk.algorithms.annotation.spec import load_spec

    doc = (
        "name: t\n"
        "layer1:\n"
        "  - axis: ptm_density\n"
        "    source: uniprot-ptm:sites\n"
        "    key: variant\n"
        "    columns: [ptm_nobs_mean, n_ptm_sites]\n"
        "    aggregate:\n"
        "      by: gene_symbol\n"
        "      to: symbol\n"
        "      count_name: n_ptm_sites\n"
        "      scores:\n"
        "        ptm_nobs: {column: n_observations, stats: [mean]}\n"
    )
    p = tmp_path / "s.yaml"
    p.write_text(doc)
    assert load_spec(p).entry("ptm_density").aggregate.count_name == "n_ptm_sites"


def test_variant_key_without_aggregate_is_rejected(tmp_path):
    from hvantk.algorithms.annotation.spec import load_spec

    doc = (
        "name: t\n"
        "layer1:\n"
        "  - {axis: x, source: s:d, key: variant, columns: [c1]}\n"
    )
    p = tmp_path / "s.yaml"
    p.write_text(doc)
    with pytest.raises(
        Exception
    ):  # jsonschema ValidationError (if/then) or ValueError (load_spec check)
        load_spec(p)


def test_aggregate_on_a_non_variant_key_is_allowed_by_schema_but_ignored(tmp_path):
    # A gene_id entry with no aggregate still parses; aggregate defaults to None.
    from hvantk.algorithms.annotation.spec import load_spec

    doc = (
        "name: t\n"
        "layer1:\n"
        "  - {axis: c, source: gnomad-metrics:metrics, key: gene_id, columns: [mis_z]}\n"
    )
    p = tmp_path / "s.yaml"
    p.write_text(doc)
    assert load_spec(p).entry("c").aggregate is None


def test_matrix_entry_parses(tmp_path):
    from hvantk.algorithms.annotation.spec import load_spec

    doc = (
        "name: t\n"
        "layer1:\n"
        "  - axis: expr\n"
        "    source: ucsc-cellbrowser:asp_2019\n"
        "    key: symbol\n"
        "    columns: [asp_cm_spec]\n"
        "    matrix:\n"
        "      group_axis: celltype\n"
        "      atlas: asp\n"
        "      tissue_tag: cardiac\n"
        "      specificity: {method: ewce_fraction, targets: [Ventricular cardiomyocytes], combine: max, name: cm_spec}\n"
    )
    p = tmp_path / "s.yaml"
    p.write_text(doc)
    e = load_spec(p).entry("expr")
    assert e.matrix.group_axis == "celltype" and e.matrix.atlas == "asp"
    assert e.matrix.specificity.method == "ewce_fraction"
    assert e.matrix.specificity.targets == ("Ventricular cardiomyocytes",)
    assert e.matrix.specificity.combine == "max"


def test_matrix_requires_symbol_key(tmp_path):
    from hvantk.algorithms.annotation.spec import load_spec

    doc = (
        "name: t\n"
        "layer1:\n"
        "  - axis: expr\n"
        "    source: s:a\n"
        "    key: gene_id\n"
        "    columns: [asp_cm_spec]\n"
        "    matrix: {group_axis: celltype, atlas: asp}\n"
    )
    p = tmp_path / "s.yaml"
    p.write_text(doc)
    with pytest.raises(ValueError, match="matrix"):
        load_spec(p)
