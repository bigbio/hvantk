import os
from pathlib import Path

import pytest
import json
from click.testing import CliRunner
from hvantk.commands import expression_atlas_downloader
from hvantk.tests.test_create_gene_annotation_tables import TEST_DIR

TEST_DIR = Path(__file__).parent.parent

@pytest.fixture
def download_path(tmpdir):
    return str(tmpdir.mkdir("downloads"))



def test_download_experiments(download_path):
    # Create a CliRunner instance
    runner = CliRunner()

    # Call the download_experiments function with the config file and download path
    config_path = TEST_DIR / "resources/expression_atlas.json"
    result = runner.invoke(expression_atlas_downloader.download_experiments, ["--config_path", config_path, "--download_path", download_path])

    # Assert that the command was successful
    assert result.exit_code == 0

    # Assert that the file was created
    file_path = os.path.join(download_path, "E-GTEX-8.condensed-sdrf.tsv")
    assert os.path.exists(file_path)

    # Assert that the file content is not empty
    assert os.path.getsize(file_path) > 0