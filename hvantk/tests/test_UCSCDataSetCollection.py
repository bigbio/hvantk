import json

from hvantk.datasets.ucsc_cell_datasets import UCSCDataSetCollection, UCSCDataset
from hvantk.utils.constants import UCSC_JSON_FILE_PATH

# Correctly create UCSCDataset objects from JSON data
def test_from_json_creates_dataset_objects(tmp_path):
    # Create a temporary JSON file with multiple datasets
    json_path = tmp_path / "test_datasets.json"
    test_data = {
        "shortLabel": "Test Label",
        "abstract": "Test Abstract",
        "inDir": "/test/dir",
        "name": "test_collection",
        "datasets": [
            {
                "shortLabel": "DS1",
                "name": "dataset1",
                "md5": "abc123",
                "tags": ["tag1", "tag2"],
                "hasFiles": ["file1.txt"],
                "body_parts": ["brain"],
                "diseases": ["none"],
                "organisms": ["human"],
                "projects": ["test_project"],
                "sampleCount": 10,
            },
            {
                "shortLabel": "DS2",
                "name": "dataset2",
                "md5": "def456",
                "tags": ["tag3"],
                "hasFiles": ["file2.txt"],
                "body_parts": ["liver"],
                "diseases": ["cancer"],
                "organisms": ["mouse"],
                "projects": ["test_project"],
                "sampleCount": 20,
            },
        ],
    }

    with open(json_path, "w") as f:
        json.dump(test_data, f)

    # Load the collection from the JSON file
    collection = UCSCDataSetCollection.from_json(str(json_path))

    # Verify the datasets were created correctly
    assert len(collection.datasets) == 2
    assert isinstance(collection.datasets[0], UCSCDataset)
    assert collection.datasets[0].name == "dataset1"
    assert collection.datasets[0].shortLabel == "DS1"
    assert collection.datasets[0].sampleCount == 10
    assert collection.datasets[1].name == "dataset2"
    assert collection.datasets[1].shortLabel == "DS2"
    assert collection.datasets[1].sampleCount == 20


# Test create UCSCDataSetCollection object from existing JSON file
def test_from_json_existing_file():
    collection = UCSCDataSetCollection.from_json(
       UCSC_JSON_FILE_PATH
    )
    assert len(collection.datasets) == 267
    assert isinstance(collection.datasets[0], UCSCDataset)
    assert collection.datasets[0].name == "cortex-dev"

    # Print the summary of the first dataset
    print(collection.datasets[0].summary())


def test_list_dataset_names():
    collection = UCSCDataSetCollection.from_json(
        UCSC_JSON_FILE_PATH
    )
    dataset_names = collection.list_dataset_names()
    print(dataset_names)
    assert isinstance(dataset_names, list)
    assert len(dataset_names) == 267
    assert "dev-brain-regions" in dataset_names
