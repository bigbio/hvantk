"""CLI wiring for `hvantk cohort`.

The parse/validate paths are pure Python and run in the fast suite; only the end-to-end
attach needs Hail.
"""
import json

import pytest
from click.testing import CliRunner

from hvantk.tools.cohort.cohort_cli import cohort_group

MANIFEST = """\
name: demo
key: gene_id
table: {table}
prior:
  column: minp
  direction: lower_is_better
cohort_axes:
  - axis: burden
    columns: [n_case_var]
"""

MANIFEST_WITH_KEY_COLUMN = """\
name: demo
key: symbol
key_column: gene
table: {table}
prior:
  column: minp
  direction: lower_is_better
"""


def _write_manifest(tmp_path, table_path):
    p = tmp_path / "cohort.yaml"
    p.write_text(MANIFEST.format(table=table_path))
    return str(p)


def _write_cohort_tsv(tmp_path):
    p = tmp_path / "cohort.tsv"
    p.write_text("gene_id\tminp\tn_case_var\nENSG1\t0.01\t5\nENSG2\t0.2\t3\n")
    return str(p)


def test_cohort_group_exposes_validate_and_attach():
    assert "validate" in cohort_group.commands
    assert "attach" in cohort_group.commands


def test_validate_reports_the_declared_columns(tmp_path):
    manifest = _write_manifest(tmp_path, _write_cohort_tsv(tmp_path))
    result = CliRunner().invoke(cohort_group, ["validate", "--cohort", manifest])
    assert result.exit_code == 0, result.output
    assert "minp" in result.output
    assert "n_case_var" in result.output


def test_validate_fails_loud_when_a_declared_column_is_missing(tmp_path):
    tsv = tmp_path / "cohort.tsv"
    tsv.write_text("gene_id\tminp\nENSG1\t0.01\n")
    manifest = _write_manifest(tmp_path, str(tsv))

    result = CliRunner().invoke(cohort_group, ["validate", "--cohort", manifest])
    assert result.exit_code != 0
    assert "n_case_var" in result.output


def test_validate_fails_loud_when_the_key_column_is_missing(tmp_path):
    """The exact real-world flaw: key='symbol' but the table's identifier column is
    named 'gene', not 'symbol'. Declaring key_column='gene' names it correctly, but
    the table below has neither -- validate must catch this instead of letting a
    later hl.import_table blow up with a Hail-internal error."""
    tsv = tmp_path / "cohort.tsv"
    tsv.write_text("symbol\tminp\nA\t0.01\n")
    p = tmp_path / "cohort.yaml"
    p.write_text(MANIFEST_WITH_KEY_COLUMN.format(table=str(tsv)))

    result = CliRunner().invoke(cohort_group, ["validate", "--cohort", str(p)])
    assert result.exit_code != 0
    assert "gene" in result.output


def test_validate_fails_loud_when_the_table_is_missing(tmp_path):
    manifest = _write_manifest(tmp_path, str(tmp_path / "nope.tsv"))
    result = CliRunner().invoke(cohort_group, ["validate", "--cohort", manifest])
    assert result.exit_code != 0
    assert "nope.tsv" in result.output


def test_validate_surfaces_a_bad_manifest_as_a_clean_error(tmp_path):
    """A schema-invalid manifest must exit non-zero with the actionable text, not a
    raw traceback. `not isinstance(..., KeyError)` alone is tautological here: schema
    validation always rejects a bad document before any `doc["table"]` lookup runs, so
    KeyError can never be raised regardless of whether `_load_manifest`'s try/except
    does anything at all -- confirmed by temporarily deleting that try/except, which
    makes THIS rewritten assertion fail (a raw jsonschema.ValidationError escapes and
    click reports a non-clean exit) while the old assertion kept passing."""
    p = tmp_path / "bad.yaml"
    p.write_text("name: demo\nkey: symbol\n")  # no table, no prior
    result = CliRunner().invoke(cohort_group, ["validate", "--cohort", str(p)])
    assert result.exit_code != 0
    assert "invalid cohort manifest" in result.output
    assert not isinstance(result.exception, KeyError), result.exception


def test_validate_succeeds_on_a_gzipped_tsv_table(tmp_path):
    """finding 3: `Path("x.tsv.gz").suffix == ".gz"`, so the old suffix-whitelist
    delimiter heuristic picked ',' and then opened the gzip bytes as UTF-8 text --
    exit 1, empty output, UnicodeDecodeError on byte 0x8b. Hail's own import_table
    already reads gz/bgz natively, so the header reader was the only blocker."""
    import gzip

    table = tmp_path / "cohort.tsv.gz"
    with gzip.open(table, "wt") as fh:
        fh.write("gene_id\tminp\tn_case_var\nENSG1\t0.01\t5\nENSG2\t0.2\t3\n")
    manifest = _write_manifest(tmp_path, str(table))

    result = CliRunner().invoke(cohort_group, ["validate", "--cohort", manifest])
    assert result.exit_code == 0, result.output
    assert "minp" in result.output
    assert "n_case_var" in result.output


@pytest.mark.hail
def test_attach_writes_a_gene_id_keyed_table_and_report(hail_session, tmp_path):
    import hail as hl

    layer1_path = str(tmp_path / "layer1.ht")
    hl.Table.parallelize(
        [
            {"gene_id": "ENSG1", "mis_z": 1.0},
            {"gene_id": "ENSG2", "mis_z": 2.0},
            {"gene_id": "ENSG3", "mis_z": 3.0},
        ],
        hl.tstruct(gene_id=hl.tstr, mis_z=hl.tfloat64),
        key=["gene_id"],
    ).write(layer1_path, overwrite=True)

    manifest = _write_manifest(tmp_path, _write_cohort_tsv(tmp_path))
    out = str(tmp_path / "attached.ht")
    report_path = str(tmp_path / "report.json")

    result = CliRunner().invoke(
        cohort_group,
        [
            "attach",
            "--cohort",
            manifest,
            "--layer1",
            layer1_path,
            "--output",
            out,
            "--report",
            report_path,
        ],
    )
    assert result.exit_code == 0, result.output

    ht = hl.read_table(out)
    assert list(ht.key) == ["gene_id"]
    assert ht.count() == 3
    assert "cohort_tested" in set(ht.row)

    report = json.loads((tmp_path / "report.json").read_text())
    assert report["n_genes"] == 3
    assert report["n_tested"] == 2
