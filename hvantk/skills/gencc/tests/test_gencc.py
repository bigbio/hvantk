"""
GenCC-specific tests.

Tests GenCC dataset class (no Hail) and GenCC streamer submitter methods
(requires Hail). Shared base-class logic is already exercised by
test_clingen_streamer.py through ClinGenGeneDiseaseTableStreamer.

The GenCC streamer itself (``hvantk.skills.gencc.streamers``) lives inside
this plugin folder as part of the ``GeneDiseaseTableStreamer`` family;
streamer tests here depend on the builder to materialize the input HT.
"""

from pathlib import Path

import pytest

# --- Dataset class tests (no Hail) ---


class TestGenCCSubmissionsDataset:
    def test_from_date_valid(self):
        from hvantk.skills.gencc.shared.datasets import GenCCSubmissionsDataset

        ds = GenCCSubmissionsDataset.from_date("2025-10-15")
        assert ds.version_date == "2025-10-15"
        assert ds.file_name == "gencc-submissions-2025-10-15.tsv"
        assert "thegencc.org" in ds.download_url

    def test_from_date_invalid_format(self):
        from hvantk.skills.gencc.shared.datasets import GenCCSubmissionsDataset

        with pytest.raises(ValueError, match="Invalid version_date format"):
            GenCCSubmissionsDataset.from_date("15-10-2025")

    def test_from_latest(self):
        from hvantk.skills.gencc.shared.datasets import GenCCSubmissionsDataset

        ds = GenCCSubmissionsDataset.from_latest()
        assert ds.version_date  # should be today's date
        assert ds.file_name.startswith("gencc-submissions-")

    def test_metadata(self):
        from hvantk.skills.gencc.shared.datasets import GenCCSubmissionsDataset

        ds = GenCCSubmissionsDataset.from_date("2025-10-15")
        meta = ds.get_metadata()
        assert meta["source"] == "GenCC (Gene Curation Coalition)"
        assert "12+" in meta["description"]


def test_raw_dir_resolves_to_the_latest_export_and_refuses_to_guess(tmp_path):
    """A second day's download leaves two dated exports in the raw directory: the
    builder must pick the newer one, and raise rather than guess when a name's date
    is not a real date or a dated file is not the same export."""
    from hvantk.core.utils.raw_inputs import resolve_raw_file
    from hvantk.skills.gencc.builder import _RAW_FILE_PATTERNS
    from hvantk.skills.gencc.shared.datasets import GenCCSubmissionsDataset

    def resolve(case, *names):
        raw_dir = tmp_path / case
        raw_dir.mkdir()
        for name in names:
            (raw_dir / name).write_text("x")
        return resolve_raw_file(
            str(raw_dir), patterns=_RAW_FILE_PATTERNS, latest_dated=True
        )

    older = GenCCSubmissionsDataset.from_date("2026-10-13").file_name
    newer = GenCCSubmissionsDataset.from_date("2026-10-14").file_name
    assert resolve("two_days", older, newer) == str(tmp_path / "two_days" / newer)

    for case, other in [
        ("invalid_date", "gencc-submissions-2026-13-45.tsv"),
        ("other_export", "gencc-other-report-2026-10-14.tsv"),
    ]:
        with pytest.raises(ValueError, match="expected exactly one file matching"):
            resolve(case, newer, other)


# --- Streamer tests (require Hail) ---

TEST_DIR = Path(__file__).parent / "testdata"


@pytest.fixture
def gencc_table_path(tmp_path):
    """Build a GenCC Hail Table from test fixture via the builder."""
    from hvantk.core.models.build_context import BuildContext
    from hvantk.skills.gencc.builder import build_gencc_submissions

    input_path = TEST_DIR / "raw/gencc/gencc_test_sample.tsv"
    output_path = tmp_path / "gencc_test.ht"
    ctx = BuildContext(
        plugin="gencc",
        dataset="gencc:submissions",
        plugin_version="test",
        source_fingerprint="sha256:test",
        builder_commit=None,
    )
    artifact = build_gencc_submissions(parsed_input=str(input_path), ctx=ctx)
    artifact.save(str(output_path))
    return str(output_path)


@pytest.mark.hail
@pytest.mark.slow
def test_gencc_submitter_summary(gencc_table_path):
    from hvantk.skills.gencc.streamers import GenCCGeneDiseaseTableStreamer

    streamer = GenCCGeneDiseaseTableStreamer(gencc_table_path, init_hail=False)
    summary = streamer.submitter_summary()
    assert "ClinGen" in summary["submitter"].values
    assert len(summary) >= 3  # ClinGen, PanelApp, Orphanet, G2P


@pytest.mark.hail
@pytest.mark.slow
def test_gencc_get_geneset_per_submitter(gencc_table_path):
    from hvantk.skills.gencc.streamers import GenCCGeneDiseaseTableStreamer

    streamer = GenCCGeneDiseaseTableStreamer(gencc_table_path, init_hail=False)
    result = streamer.get_geneset_per_submitter()
    assert "ClinGen" in result
    assert "BRCA1" in result["ClinGen"]


@pytest.mark.hail
@pytest.mark.slow
def test_gencc_consensus_genes(gencc_table_path):
    from hvantk.skills.gencc.streamers import GenCCGeneDiseaseTableStreamer

    streamer = GenCCGeneDiseaseTableStreamer(gencc_table_path, init_hail=False)
    # BRCA1, BRCA2, TP53, PTEN have multiple submitters in test data
    consensus = streamer.consensus_genes(min_submitters=2)
    assert "BRCA1" in consensus
    assert "BRCA2" in consensus
