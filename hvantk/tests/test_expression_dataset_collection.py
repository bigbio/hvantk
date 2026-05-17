import json

from hvantk.skills.expression_atlas.shared.datasets import (
    ExpressionAtlasDatasetCollection,
    ExpressionAtlasDataset,
)
from hvantk.core.constants import EXPRESSION_ATLAS_JSON_FILE_PATH


def test_from_json_creates_dataset_objects(tmp_path):
    """Test that the from_json method correctly creates dataset objects from JSON data."""
    # Create a temporary JSON file with multiple datasets
    json_path = tmp_path / "test_datasets.json"
    test_data = [
        {
            "title": "Dataset 1",
            "accession": "E-MTAB-1234",
            "type": "RNA-Seq mRNA baseline",
            "pubmedid": "12345678",
            "description": "Test description 1",
            "files": [
                {"type": "tpm", "name": "dataset1.tpm.txt"},
                {"type": "sdrf", "name": "dataset1.sdrf.txt"},
            ],
        },
        {
            "title": "Dataset 2",
            "accession": "E-MTAB-5678",
            "type": "RNA-Seq differential",
            "description": "Test description 2",
            "files": [
                {"type": "fpkm", "name": "dataset2.fpkm.txt"},
                {"type": "sdrf", "name": "dataset2.sdrf.txt"},
            ],
        },
    ]

    with open(json_path, "w") as f:
        json.dump(test_data, f)

    # Load the collection from the JSON file
    collection = ExpressionAtlasDatasetCollection.from_json(str(json_path))

    # Verify the datasets were created correctly
    assert len(collection.datasets) == 2
    assert isinstance(collection.datasets[0], ExpressionAtlasDataset)
    assert collection.datasets[0].accession == "E-MTAB-1234"
    assert collection.datasets[0].title == "Dataset 1"
    assert collection.datasets[0].type == "RNA-Seq mRNA baseline"
    assert collection.datasets[0].pubmedid == "12345678"
    assert len(collection.datasets[0].files) == 2
    assert collection.datasets[1].accession == "E-MTAB-5678"
    assert collection.datasets[1].title == "Dataset 2"
    assert collection.datasets[1].type == "RNA-Seq differential"
    assert collection.datasets[1].pubmedid is None


def test_from_json_existing_file():
    """Test that from_json correctly loads from an existing file."""
    collection = ExpressionAtlasDatasetCollection.from_json(
        EXPRESSION_ATLAS_JSON_FILE_PATH
    )
    assert len(collection.datasets) > 0
    assert isinstance(collection.datasets[0], ExpressionAtlasDataset)

    # Print the summary of the first dataset
    print(collection.datasets[1].summary())


def test_get_dataset_by_accession():
    """Test the get_by_accession method."""
    # Create a collection with test datasets
    dataset1 = ExpressionAtlasDataset(
        title="Dataset 1", accession="E-MTAB-1234", type="RNA-Seq mRNA baseline"
    )
    dataset2 = ExpressionAtlasDataset(
        title="Dataset 2", accession="E-MTAB-5678", type="RNA-Seq differential"
    )
    collection = ExpressionAtlasDatasetCollection(datasets=[dataset1, dataset2])

    # Test finding an existing dataset
    found_dataset = collection.get_by_accession("E-MTAB-1234")
    assert found_dataset is not None
    assert found_dataset.title == "Dataset 1"

    # Test finding a non-existent dataset
    not_found_dataset = collection.get_by_accession("non-existent")
    assert not_found_dataset is None


def test_list_dataset_accessions():
    """Test the list_dataset_accessions method."""
    collection = ExpressionAtlasDatasetCollection.from_json(
        EXPRESSION_ATLAS_JSON_FILE_PATH
    )
    accessions = collection.list_dataset_accessions()
    assert len(accessions) > 0
    assert "E-GTEX-8" in accessions
    assert "E-MTAB-6782" in accessions
