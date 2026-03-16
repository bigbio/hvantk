import json
from unittest.mock import patch, MagicMock

from hvantk.datasets.ucsc_cell_datasets import UCSCDataSetCollection, UCSCDataset
from hvantk.core.constants import UCSC_JSON_FILE_PATH


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
    collection = UCSCDataSetCollection.from_json(UCSC_JSON_FILE_PATH)
    assert len(collection.datasets) == 267
    assert isinstance(collection.datasets[0], UCSCDataset)
    assert collection.datasets[0].name == "cortex-dev"

    # Print the summary of the first dataset
    print(collection.datasets[0].summary())


def test_list_dataset_names():
    collection = UCSCDataSetCollection.from_json(UCSC_JSON_FILE_PATH)
    dataset_names = collection.list_dataset_names()
    print(dataset_names)
    assert isinstance(dataset_names, list)
    assert len(dataset_names) == 267
    assert "dev-brain-regions" in dataset_names


# --- Phase 2.2: search() tests ---


class TestSearch:
    """Tests for UCSCDataSetCollection.search()."""

    def test_search_by_name(self):
        """Search matches dataset name."""
        collection = UCSCDataSetCollection.from_json(UCSC_JSON_FILE_PATH)
        results = collection.search("adultPancreas")
        assert len(results.datasets) >= 1
        assert any(ds.name == "adultPancreas" for ds in results.datasets)

    def test_search_case_insensitive(self):
        """Search is case-insensitive."""
        collection = UCSCDataSetCollection.from_json(UCSC_JSON_FILE_PATH)
        results_lower = collection.search("pancreas")
        results_upper = collection.search("PANCREAS")
        assert len(results_lower.datasets) == len(results_upper.datasets)

    def test_search_by_body_part(self):
        """Search matches body_parts facet."""
        collection = UCSCDataSetCollection.from_json(UCSC_JSON_FILE_PATH)
        results = collection.search("heart")
        assert len(results.datasets) >= 1

    def test_search_no_results(self):
        """Search with no matches returns empty collection."""
        collection = UCSCDataSetCollection.from_json(UCSC_JSON_FILE_PATH)
        results = collection.search("zzz_nonexistent_zzz")
        assert len(results.datasets) == 0

    def test_search_returns_collection_type(self):
        """Search returns a UCSCDataSetCollection."""
        collection = UCSCDataSetCollection.from_json(UCSC_JSON_FILE_PATH)
        results = collection.search("brain")
        assert isinstance(results, UCSCDataSetCollection)


# --- Phase 1.1: fetch_children() tests ---


class TestFetchChildren:
    """Tests for UCSCDataset.fetch_children()."""

    def test_fetch_children_not_collection(self):
        """Leaf datasets return empty list without network call."""
        ds = UCSCDataset(shortLabel="Leaf", name="leaf-ds", md5="", isCollection=False)
        assert ds.fetch_children() == []

    def test_fetch_children_caches_result(self):
        """Second call returns cached children, no extra HTTP request."""
        ds = UCSCDataset(shortLabel="Coll", name="test-coll", md5="", isCollection=True)
        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "datasets": [
                {"shortLabel": "Child1", "name": "test-coll/child1", "sampleCount": 100}
            ]
        }
        with patch("requests.get", return_value=mock_resp) as mock_get:
            children1 = ds.fetch_children()
            children2 = ds.fetch_children()
            assert mock_get.call_count == 1
        assert len(children1) == 1
        assert children1 is children2

    def test_fetch_children_parses_response(self):
        """Children are parsed into UCSCDataset objects."""
        ds = UCSCDataset(
            shortLabel="HOC", name="hoc", md5="", isCollection=True, datasetCount=2
        )
        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "datasets": [
                {
                    "shortLabel": "All Heart",
                    "name": "hoc/all-heart",
                    "sampleCount": 142946,
                    "organisms": ["Human"],
                    "body_parts": ["heart"],
                },
                {
                    "shortLabel": "Blood",
                    "name": "hoc/blood",
                    "sampleCount": 12345,
                },
            ]
        }
        with patch("requests.get", return_value=mock_resp):
            children = ds.fetch_children()
        assert len(children) == 2
        assert children[0].name == "hoc/all-heart"
        assert children[0].sampleCount == 142946
        assert children[1].name == "hoc/blood"
        assert isinstance(children[0], UCSCDataset)

    def test_fetch_children_network_failure(self):
        """Network failure returns empty list, no exception raised."""
        ds = UCSCDataset(shortLabel="Coll", name="bad-coll", md5="", isCollection=True)
        with patch("requests.get", side_effect=ConnectionError("timeout")):
            children = ds.fetch_children()
        assert children == []
        assert ds.children is None  # Not cached on failure
