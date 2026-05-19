"""
GenCC-specific tests.

Tests GenCC dataset class (no Hail) and GenCC streamer submitter methods
(requires Hail). Shared base-class logic is already exercised by
test_clingen_streamer.py through ClinGenStreamer.

The GenCC streamer itself (``hvantk.data.gencc_streamer``) intentionally
lives outside this plugin folder because it is part of the cross-cutting
``GeneDiseaseValidityStreamer`` family; the streamer tests here live with
the builder because they depend on it to materialize the input HT.
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


# --- Streamer tests (require Hail) ---

TEST_DIR = Path(__file__).parent / "testdata"


@pytest.fixture
def gencc_table_path(tmp_path):
    """Build a GenCC Hail Table from test fixture (gene_disease_submitter keying)."""
    from hvantk.skills.gencc.builder import create_gencc_submissions_tb

    input_path = TEST_DIR / "raw/gencc/gencc_test_sample.tsv"
    output_path = tmp_path / "gencc_test.ht"
    create_gencc_submissions_tb(
        input_path=str(input_path),
        output_path=str(output_path),
        key_by="gene_disease_submitter",
        overwrite=True,
    )
    return str(output_path)


@pytest.mark.hail
@pytest.mark.slow
def test_gencc_submitter_summary(gencc_table_path):
    from hvantk.core.streamers.gencc import GenCCStreamer

    streamer = GenCCStreamer(gencc_table_path)
    summary = streamer.submitter_summary()
    assert "ClinGen" in summary["submitter"].values
    assert len(summary) >= 3  # ClinGen, PanelApp, Orphanet, G2P


@pytest.mark.hail
@pytest.mark.slow
def test_gencc_get_geneset_per_submitter(gencc_table_path):
    from hvantk.core.streamers.gencc import GenCCStreamer

    streamer = GenCCStreamer(gencc_table_path)
    result = streamer.get_geneset_per_submitter()
    assert "ClinGen" in result
    assert "BRCA1" in result["ClinGen"]


@pytest.mark.hail
@pytest.mark.slow
def test_gencc_consensus_genes(gencc_table_path):
    from hvantk.core.streamers.gencc import GenCCStreamer

    streamer = GenCCStreamer(gencc_table_path)
    # BRCA1, BRCA2, TP53, PTEN have multiple submitters in test data
    consensus = streamer.consensus_genes(min_submitters=2)
    assert "BRCA1" in consensus
    assert "BRCA2" in consensus
