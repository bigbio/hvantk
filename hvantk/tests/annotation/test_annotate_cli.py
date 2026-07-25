"""CLI tests for ``hvantk annotate compose``.

The malformed-``--prepared`` and axis-validation checks are pure Python (no spec,
no Hail) and run in the fast suite; the end-to-end invocation writes/reads real Hail
tables and is marked ``@pytest.mark.hail``.
"""
from __future__ import annotations

import json

import click
import pytest
from click.testing import CliRunner

from hvantk.tools.annotation.annotate_cli import (
    _parse_prepared_option,
    _validate_prepared_axes,
    annotate_group,
)


def test_annotate_group_exposes_compose():
    assert "compose" in annotate_group.commands


# --------------------------------------------------------------------------
# _parse_prepared_option -- pure Python, no Hail
# --------------------------------------------------------------------------


def test_parse_prepared_option_builds_axis_to_path_dict():
    parsed = _parse_prepared_option(("constraint=/tmp/c.ht", "gevir=/tmp/g.ht"))
    assert parsed == {"constraint": "/tmp/c.ht", "gevir": "/tmp/g.ht"}


def test_parse_prepared_option_rejects_a_value_with_no_equals_sign():
    with pytest.raises(click.UsageError, match="malformed"):
        _parse_prepared_option(("constraint-only-path",))


def test_parse_prepared_option_rejects_an_empty_axis_or_path():
    with pytest.raises(click.UsageError, match="malformed"):
        _parse_prepared_option(("=/tmp/c.ht",))
    with pytest.raises(click.UsageError, match="malformed"):
        _parse_prepared_option(("constraint=",))


def test_parse_prepared_option_rejects_a_duplicate_axis():
    # A repeated axis must not silently last-win -- a typo'd repeat would otherwise
    # quietly compose the wrong table.
    with pytest.raises(click.UsageError, match="twice"):
        _parse_prepared_option(("constraint=/tmp/a.ht", "constraint=/tmp/b.ht"))


# --------------------------------------------------------------------------
# _validate_prepared_axes -- pure Python, no Hail
# --------------------------------------------------------------------------


class _FakeEntry:
    def __init__(self, axis):
        self.axis = axis


class _FakeSpec:
    def __init__(self, name, axes):
        self.name = name
        self.layer1 = tuple(_FakeEntry(a) for a in axes)


def test_validate_prepared_axes_passes_when_they_match_exactly():
    spec = _FakeSpec("t", ["constraint", "gevir"])
    _validate_prepared_axes(
        spec, {"constraint": "/tmp/c.ht", "gevir": "/tmp/g.ht"}
    )  # must not raise


def test_validate_prepared_axes_raises_on_a_missing_spec_axis():
    spec = _FakeSpec("t", ["constraint", "gevir"])
    with pytest.raises(click.UsageError, match="gevir"):
        _validate_prepared_axes(spec, {"constraint": "/tmp/c.ht"})


def test_validate_prepared_axes_raises_on_an_unknown_axis():
    spec = _FakeSpec("t", ["constraint"])
    with pytest.raises(click.UsageError, match="unknown_axis"):
        _validate_prepared_axes(
            spec, {"constraint": "/tmp/c.ht", "unknown_axis": "/tmp/u.ht"}
        )


# --------------------------------------------------------------------------
# End-to-end CLI invocation -- writes/reads real Hail tables.
# --------------------------------------------------------------------------


def _write_spine(path):
    import hail as hl

    ht = hl.Table.parallelize(
        [
            {"gene_id": "ENSG1", "gene_name": "A"},
            {"gene_id": "ENSG2", "gene_name": "B"},
            {"gene_id": "ENSG3", "gene_name": "C"},
        ],
        hl.tstruct(gene_id=hl.tstr, gene_name=hl.tstr),
        key=["gene_id"],
    )
    ht.write(path, overwrite=True)


def _write_constraint_axis(path):
    import hail as hl

    # ENSG1/ENSG2 only -- ENSG3 absent, exercising the no-imputation contract.
    ht = hl.Table.parallelize(
        [
            {"gene_id": "ENSG1", "mis_z": 2.0},
            {"gene_id": "ENSG2", "mis_z": -1.0},
        ],
        hl.tstruct(gene_id=hl.tstr, mis_z=hl.tfloat64),
        key=["gene_id"],
    )
    ht.write(path, overwrite=True)


def _write_gevir_axis(path):
    import hail as hl

    ht = hl.Table.parallelize(
        [
            {"gene_id": "ENSG2", "gevir_pct": 0.5},
            {"gene_id": "ENSG3", "gevir_pct": 0.9},
        ],
        hl.tstruct(gene_id=hl.tstr, gevir_pct=hl.tfloat64),
        key=["gene_id"],
    )
    ht.write(path, overwrite=True)


def _write_spec(path):
    path.write_text(
        "name: test-spec\n"
        "layer1:\n"
        "  - {axis: constraint, source: gnomad-metrics:metrics, key: gene_id, "
        "columns: [mis_z]}\n"
        "  - {axis: gevir, source: gevir:table, key: gene_id, columns: [gevir_pct]}\n"
    )


@pytest.mark.hail
def test_compose_cmd_writes_a_gene_id_keyed_matrix_and_manifest(hail_session, tmp_path):
    import hail as hl

    spec_path = tmp_path / "spec.yaml"
    _write_spec(spec_path)

    spine_path = str(tmp_path / "spine.ht")
    _write_spine(spine_path)

    constraint_path = str(tmp_path / "constraint.ht")
    _write_constraint_axis(constraint_path)

    gevir_path = str(tmp_path / "gevir.ht")
    _write_gevir_axis(gevir_path)

    output_path = str(tmp_path / "matrix.ht")
    manifest_path = str(tmp_path / "manifest.json")

    runner = CliRunner()
    result = runner.invoke(
        annotate_group,
        [
            "compose",
            "--spec",
            str(spec_path),
            "--spine",
            spine_path,
            "--prepared",
            f"constraint={constraint_path}",
            "--prepared",
            f"gevir={gevir_path}",
            "--output",
            output_path,
            "--manifest",
            manifest_path,
        ],
    )

    assert result.exit_code == 0, result.output

    ht = hl.read_table(output_path)
    assert list(ht.key) == ["gene_id"]
    assert ht.count() == 3
    assert set(ht.row) == {
        "gene_id",
        "gene_name",
        "mis_z",
        "gevir_pct",
        "constraint_present",
        "gevir_present",
    }

    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["n_genes"] == 3
    assert set(manifest["axes"]) == {"constraint", "gevir"}


@pytest.mark.hail
def test_compose_cmd_defaults_manifest_path_alongside_output(hail_session, tmp_path):
    spec_path = tmp_path / "spec.yaml"
    _write_spec(spec_path)

    spine_path = str(tmp_path / "spine.ht")
    _write_spine(spine_path)

    constraint_path = str(tmp_path / "constraint.ht")
    _write_constraint_axis(constraint_path)

    gevir_path = str(tmp_path / "gevir.ht")
    _write_gevir_axis(gevir_path)

    output_path = str(tmp_path / "matrix.ht")

    runner = CliRunner()
    result = runner.invoke(
        annotate_group,
        [
            "compose",
            "--spec",
            str(spec_path),
            "--spine",
            spine_path,
            "--prepared",
            f"constraint={constraint_path}",
            "--prepared",
            f"gevir={gevir_path}",
            "--output",
            output_path,
        ],
    )

    assert result.exit_code == 0, result.output
    assert (tmp_path / "matrix.ht.manifest.json").exists()


def test_compose_cmd_errors_clearly_on_a_missing_axis(tmp_path):
    spec_path = tmp_path / "spec.yaml"
    _write_spec(spec_path)

    runner = CliRunner()
    result = runner.invoke(
        annotate_group,
        [
            "compose",
            "--spec",
            str(spec_path),
            "--spine",
            str(tmp_path / "spine.ht"),
            "--prepared",
            f"constraint={tmp_path / 'constraint.ht'}",
            "--output",
            str(tmp_path / "matrix.ht"),
        ],
    )

    assert result.exit_code != 0
    assert not isinstance(result.exception, KeyError)
    assert "gevir" in result.output


def test_compose_cmd_errors_clearly_on_an_unknown_axis(tmp_path):
    spec_path = tmp_path / "spec.yaml"
    _write_spec(spec_path)

    runner = CliRunner()
    result = runner.invoke(
        annotate_group,
        [
            "compose",
            "--spec",
            str(spec_path),
            "--spine",
            str(tmp_path / "spine.ht"),
            "--prepared",
            f"constraint={tmp_path / 'constraint.ht'}",
            "--prepared",
            f"gevir={tmp_path / 'gevir.ht'}",
            "--prepared",
            f"expression={tmp_path / 'expression.ht'}",
            "--output",
            str(tmp_path / "matrix.ht"),
        ],
    )

    assert result.exit_code != 0
    assert "expression" in result.output


def test_compose_cmd_errors_clearly_on_a_malformed_prepared_value(tmp_path):
    spec_path = tmp_path / "spec.yaml"
    _write_spec(spec_path)

    runner = CliRunner()
    result = runner.invoke(
        annotate_group,
        [
            "compose",
            "--spec",
            str(spec_path),
            "--spine",
            str(tmp_path / "spine.ht"),
            "--prepared",
            "not-a-key-value-pair",
            "--output",
            str(tmp_path / "matrix.ht"),
        ],
    )

    assert result.exit_code != 0
    assert "malformed" in result.output
