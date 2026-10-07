"""Round-trip test for the PeptideAtlas phospho plugin: parse + build.

Runs the real ``parse_raw_dir`` over a synthetic raw-build zip fixture, then
``build_peptideatlas_phospho`` over the resulting intermediate TSV, and
asserts the schema and a sample of rows against committed snapshots -- so a
change in either the raw-zip parsing contract (table joins, phospho offset
extraction, canonical / DECOY_ / CONTAM_ filtering, observation-count
aggregation) or the intermediate-TSV -> ``AnnotationTable`` contract shows up
as an explicit diff rather than silently.

Fixture provenance -- READ THIS BEFORE TRUSTING THE NUMBERS
-------------------------------------------------------------
``testdata/raw/peptideatlas-phospho/atlas_build_606-synthetic.tsv.zip`` is a
**synthetic miniature raw build**, not a truncation of a real PeptideAtlas
build: the real ``atlas_build_*.tsv.zip`` is ~549 MB (``content_length`` in
the committed ``tests/drift_fingerprint.json``) and is not (and should not
be) vendored into this repo. Its five tables carry the REAL column headers
(verbatim, same order) captured from a real human phospho build (202512 /
606) on 2026-10-06, and the REAL modification notation
(``<residue>[Phospho]``), but every row is fabricated -- see
``testdata/raw/peptideatlas-phospho/README.md`` for the full row-purpose
table and the "what's real, what's fabricated" breakdown. In short:

  - ``P04637`` (fake seq/counts), canonical: an aggregated T@14 site (two
    peptides, 42+17=59 observations), a multi-site peptide giving S@32 and
    Y@35 (30 observations each), and one arm of a multi-mapping peptide
    giving S@53 (8 observations).
  - ``P99999``, canonical: the other arm of that same multi-mapping peptide,
    S@23 (8 observations).
  - ``Q99999``, a non-canonical isoform (``presence_level_id`` = 3): must be
    filtered out by the canonical check.
  - ``DECOY_FAKE1`` / ``CONTAM_FAKE1``: must be filtered out by accession
    prefix.
  - One peptide carries only a non-phospho modification
    (``[Carbamidomethyl]``): must contribute zero output rows.
  - One peptide (``peptide_instance`` 209, seen 100 times) in two phospho
    forms (12 and 8 observations) plus its unmodified form (80): its site
    P04637 S@42 counts 12 + 8 = 20, each form's own observations (#425).

Regenerate the zip itself only by re-deriving it from the row design above
(and in the README) -- no generator script is committed to this repo.

What this test covers
----------------------
Both stages of the plugin's own contract: ``parse_raw_dir`` (raw zip ->
intermediate TSV -- table joins, phospho-offset extraction, canonical /
DECOY_ / CONTAM_ filtering, observation-count aggregation; see
``shared/datasets.py::parse_peptideatlas_zip``) and
``build_peptideatlas_phospho`` (intermediate TSV -> ``AnnotationTable``; a
thin ``pd.read_csv`` + ``AnnotationTable.from_pandas`` wrap -- see
``builder.py``'s module docstring). ``test_phospho.py`` separately
unit-tests ``parse_peptideatlas_zip`` / ``_extract_phospho_offsets`` against
hand-written mock tables with their own notation and edge cases (numeric-mass
brackets, N-terminal labels, canonical filtering in isolation, etc.); this
test is the one committed fixture that exercises the full raw-zip ->
built-artifact path the plugin ships end to end.

Why this does not use ``hvantk.tests._snapshot_utils.regenerate_snapshots``
------------------------------------------------------------------------------
``build_peptideatlas_phospho`` returns a **pandas-backed** ``AnnotationTable``
(``peptideatlas`` has no native Hail/AnnData representation -- see
``builder.py``). ``regenerate_snapshots`` only has two dispatch branches: one
for ``anndata.AnnData`` and a Hail-table fallback that re-reads whatever was
written to an ``output_path``. Neither is a true fit: this artifact is not
AnnData, and forcing it through the Hail branch would only work via
``AnnotationTable.to_hail()`` (which drops the pandas dtypes and requires a
full Hail/Spark session, `@pytest.mark.hail`, just to validate a builder that
never touches Hail -- see ``hail_session`` cost in ``hvantk/tests/conftest.py``
and CLAUDE.md's "one contract per layer" guidance). ``AnnotationTable`` already
computes its own backend-agnostic ``.schema`` dict
(``hvantk/core/models/annotation_table.py::_normalize_dtype``), so this test
uses that directly plus a small local helper mirroring ``collect_sample_rows``'
``{"key": ..., "row": ...}`` shape, and still uses ``_snapshot_utils.load_snapshot``
for the read side. This keeps the test Hail-free, matching every other
pure-pandas plugin path and the default (non-``-m hail``) pytest selection.

Regenerate after an intentional change:
    pytest hvantk/skills/peptideatlas/phospho/tests/test_builder.py --regenerate-snapshots
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from hvantk.tests._snapshot_utils import load_snapshot

# A directory, not a file: parse_raw_dir expects exactly one
# atlas_build_*.tsv.zip inside it (see shared/datasets.py::parse_raw_dir).
FIXTURE_DIR = str(
    Path("hvantk/skills/peptideatlas/phospho/tests/testdata/raw/peptideatlas-phospho")
)
SNAPSHOT_DIR = Path("hvantk/skills/peptideatlas/phospho/tests/snapshots")

EXPECTED_COLUMNS = {
    "accession",
    "gene_symbol",
    "position",
    "description",
    "amino_acid",
    "ensembl_xrefs",
    "sequence_length",
    "n_observations",
    "source_db",
    "evidence_type",
}

# (accession, position) is the composite identity parse_peptideatlas_zip already
# dedupes/aggregates phospho sites on -- the natural row key for this table, even
# though the pandas-backed AnnotationTable itself carries no formal key concept.
KEY_FIELDS = ("accession", "position")

# Keys into the fabricated fixture (see README.md's row-purpose table beside the
# zip), chosen to cover four different code paths in one snapshot: an aggregated
# multi-peptide site, a site from a one-peptide/two-site extraction, one arm of a
# multi-mapping peptide, and a site seen in two phospho forms of one peptide
# (counted from each form's own observations, #425). Row selection below raises KeyError for any key
# absent from the table, so these cannot be invented.
SAMPLE_KEYS = [
    {"accession": "P04637", "position": "14"},
    {"accession": "P04637", "position": "35"},
    {"accession": "P99999", "position": "23"},
    {"accession": "P04637", "position": "42"},
]


def _fake_ctx():
    """Deterministic BuildContext, so provenance never perturbs a snapshot."""
    from hvantk.core.models.build_context import BuildContext

    return BuildContext(
        plugin="peptideatlas",
        dataset="peptideatlas:phospho",
        plugin_version="test",
        source_fingerprint="sha256:test",
        builder_commit=None,
    )


def _build(parsed_input: str):
    """Run the builder under test against an already-parsed intermediate TSV."""
    from hvantk.skills.peptideatlas.phospho.builder import build_peptideatlas_phospho

    return build_peptideatlas_phospho(parsed_input=parsed_input, ctx=_fake_ctx())


def _schema_to_dict(artifact, df) -> dict:
    """JSON-stable schema dict: ordered columns + AnnotationTable's normalized dtypes."""
    return {
        "n_rows": int(len(df)),
        "columns": list(df.columns),
        "dtypes": artifact.schema,
    }


def _sample_rows(df, keys: list[dict]) -> list[dict]:
    """Select rows by (accession, position), in the order of ``keys``.

    Mirrors ``_snapshot_utils.collect_sample_rows``'s ``{"key": ..., "row": ...}``
    shape and its fail-loudly-on-missing-key behavior (a silently dropped key would
    produce a confusing snapshot diff instead of an explicit error).
    """
    out = []
    for k in keys:
        mask = None
        for field in KEY_FIELDS:
            field_mask = df[field] == k[field]
            mask = field_mask if mask is None else (mask & field_mask)
        matches = df[mask]
        if matches.empty:
            available = list(zip(df["accession"], df["position"]))[:5]
            raise KeyError(
                f"requested key {k!r} not found in table; "
                f"available keys: {available} (showing up to 5)"
            )
        row = matches.iloc[0]
        # pd.read_csv(..., dtype=str) still applies pandas' default NA detection, so an
        # empty field (e.g. ensembl_xrefs, always "" in write_intermediate_tsv today)
        # comes back as float NaN inside an otherwise-str column, not "". Normalize to
        # None: NaN != NaN in Python, which would make every round trip "drift" even
        # when nothing changed, and None survives a JSON write/read unchanged.
        out.append(
            {
                "key": {f: row[f] for f in KEY_FIELDS},
                "row": {
                    c: (None if pd.isna(row[c]) else row[c])
                    for c in df.columns
                    if c not in KEY_FIELDS
                },
            }
        )
    return out


def test_peptideatlas_phospho_snapshot_round_trip(tmp_path, regenerate_snapshots):
    """Parse the raw fixture zip, build peptideatlas:phospho from it, and assert
    schema and sample-row stability against committed snapshots."""
    from hvantk.skills.peptideatlas.phospho.shared.datasets import parse_raw_dir

    intermediate_path = parse_raw_dir(FIXTURE_DIR, str(tmp_path / "intermediate.tsv"))
    artifact = _build(intermediate_path)
    df = artifact.to_pandas()

    # Structural assertions, independent of the snapshot comparison.
    assert len(df) == 6, (
        "fixture has exactly six phospho sites (see the 'expected parse "
        "output' table in README.md beside the fixture zip)"
    )
    assert set(df.columns) == EXPECTED_COLUMNS, (
        "intermediate-TSV column set drifted from the 10-column contract in "
        "shared/datasets.py::_TSV_COLUMNS"
    )

    if regenerate_snapshots:
        SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
        schema = _schema_to_dict(artifact, df)
        (SNAPSHOT_DIR / "schema.json").write_text(
            json.dumps(schema, indent=2, sort_keys=True)
        )
        rows = _sample_rows(df, SAMPLE_KEYS)
        (SNAPSHOT_DIR / "sample_rows.json").write_text(
            json.dumps(rows, indent=2, sort_keys=True)
        )
        pytest.skip(
            "Snapshots regenerated; rerun without --regenerate-snapshots to assert."
        )

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    assert _schema_to_dict(artifact, df) == expected_schema, (
        "peptideatlas:phospho schema drifted from snapshot"
    )

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    assert _sample_rows(df, SAMPLE_KEYS) == expected_rows, (
        "peptideatlas:phospho sample rows drifted from snapshot"
    )
