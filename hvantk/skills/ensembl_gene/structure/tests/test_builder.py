"""Snapshot round-trip for the ensembl-gene:structure builder.

Regenerate after an intentional change:
    pytest hvantk/skills/ensembl_gene/structure/tests/test_builder.py -m hail \
        --regenerate-snapshots
"""
from __future__ import annotations

from pathlib import Path

import pytest

from hvantk.tests._snapshot_utils import (
    collect_sample_rows,
    hail_schema_to_dict,
    load_snapshot,
    phase_b_snapshot_adapter,
)
from hvantk.tests._snapshot_utils import regenerate_snapshots as regenerate_snapshots_fn

FIXTURE = (
    "hvantk/skills/ensembl_gene/structure/tests/testdata/raw/"
    "ensembl-structure/mini.gtf"
)
SNAPSHOT_DIR = Path("hvantk/skills/ensembl_gene/structure/tests/snapshots")

SAMPLE_KEYS = [
    {"gene_id": "ENSG00000000001"},
    {"gene_id": "ENSG00000000002"},
    {"gene_id": "ENSG00000000003"},
]


@pytest.mark.hail
def test_structure_snapshot_round_trip(hail_session, tmp_path, regenerate_snapshots):
    import hail as hl

    from hvantk.skills.ensembl_gene.structure.builder import (
        build_ensembl_gene_structure,
    )

    builder = phase_b_snapshot_adapter(
        build_ensembl_gene_structure, "ensembl-gene:structure"
    )

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=builder,
            fixture_path=FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            keys=SAMPLE_KEYS,
        )
        pytest.skip("Snapshots regenerated; rerun without the flag to assert.")

    output_path = str(tmp_path / "structure.ht")
    builder(input_path=FIXTURE, output_path=output_path)
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    assert hail_schema_to_dict(ht) == expected_schema, "structure schema drifted"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    assert collect_sample_rows(ht, keys=SAMPLE_KEYS) == expected_rows


@pytest.mark.hail
def test_gene_id_is_the_key_and_is_unique(hail_session, tmp_path):
    import hail as hl

    from hvantk.skills.ensembl_gene.structure.builder import (
        build_ensembl_gene_structure,
    )

    builder = phase_b_snapshot_adapter(
        build_ensembl_gene_structure, "ensembl-gene:structure"
    )
    output_path = str(tmp_path / "structure.ht")
    builder(input_path=FIXTURE, output_path=output_path)
    ht = hl.read_table(output_path)

    assert list(ht.key) == ["gene_id"]
    assert ht.count() == ht.distinct().count()


@pytest.mark.hail
def test_protein_coding_only_param_filters(hail_session, tmp_path):
    import hail as hl

    from hvantk.skills.ensembl_gene.structure.builder import (
        build_ensembl_gene_structure,
    )

    builder = phase_b_snapshot_adapter(
        build_ensembl_gene_structure, "ensembl-gene:structure"
    )
    output_path = str(tmp_path / "structure_pc.ht")
    builder(input_path=FIXTURE, output_path=output_path, protein_coding_only=True)
    ht = hl.read_table(output_path)

    assert ht.count() == 2  # the lncRNA gene is dropped


@pytest.mark.hail
def test_builder_accepts_a_raw_directory(hail_session, tmp_path):
    """End-to-end guard for the reprocess path: `hvantk reprocess` hands the builder the
    raw directory (no parse stage), not the GTF file. The build must succeed against a
    directory containing the downloaded GTF, not crash with IsADirectoryError.

    The downloaded artifact is gzipped (``...gtf.gz``), so the fixture is gzipped into the
    committed filename here -- which also exercises the parser's gzip branch."""
    import gzip

    import hail as hl

    from hvantk.resources.ensembl_release import ENSEMBL_GTF_FILENAME
    from hvantk.skills.ensembl_gene.structure.builder import (
        build_ensembl_gene_structure,
    )

    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    with open(FIXTURE, "rb") as src, gzip.open(
        raw_dir / ENSEMBL_GTF_FILENAME, "wb"
    ) as dst:
        dst.write(src.read())  # what the downloader writes: a gzipped GTF

    builder = phase_b_snapshot_adapter(
        build_ensembl_gene_structure, "ensembl-gene:structure"
    )
    output_path = str(tmp_path / "structure_from_dir.ht")
    builder(
        input_path=str(raw_dir), output_path=output_path
    )  # reprocess passes the dir
    ht = hl.read_table(output_path)

    assert ht.count() == 3
    assert list(ht.key) == ["gene_id"]
