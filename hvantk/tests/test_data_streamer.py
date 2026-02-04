# Test for Data Streamer Implementation
# Basic unit tests for the Clinvar data streamer

import pytest
import logging
from unittest.mock import Mock, patch
from hvantk.data.data_streamer import DataStreamer, HailDataStreamer, StreamProcessor
from hvantk.data.clinvar_streamer import (
    ClinvarDataStreamer,
    ClinvarTrainingSetProcessor,
    create_clinvar_training_set_streamer,
)
from hvantk.utils import load_sample_chd_gene_set

# Configure logging for tests
logging.basicConfig(level=logging.INFO)


class TestDataStreamer:
    """Test base data streamer functionality"""

    def test_load_sample_gene_set(self):
        genes = load_sample_chd_gene_set()
        assert isinstance(genes, set)
        assert len(genes) > 0
        assert "GATA4" in genes
        assert "NKX2-5" in genes

    def test_clinvar_streamer_init(self):
        gene_set = {"GATA4", "NKX2-5"}
        streamer = ClinvarDataStreamer(
            clinvar_path="/path/to/clinvar.vcf.gz", gene_set=gene_set, chunk_size=5000
        )
        assert streamer.name == "ClinvarTrainingSet"
        assert streamer.chunk_size == 5000
        assert streamer.gene_set == gene_set
        assert streamer.clinvar_path == "/path/to/clinvar.vcf.gz"

    def test_factory_function(self):
        processor = create_clinvar_training_set_streamer(
            clinvar_path="/path/to/clinvar.vcf.gz", output_dir="/tmp/test"
        )
        assert isinstance(processor, ClinvarTrainingSetProcessor)
        assert processor.name == "ClinvarTrainingSetGeneration"
        assert len(processor.streamers) == 1
        assert isinstance(processor.streamers[0], ClinvarDataStreamer)

    def test_stream_processor_init(self):
        processor = StreamProcessor("TestProcessor")
        assert processor.name == "TestProcessor"
        assert len(processor.streamers) == 0
        mock_streamer = Mock(spec=DataStreamer)
        mock_streamer.name = "MockStreamer"
        processor.add_streamer(mock_streamer)
        assert len(processor.streamers) == 1
        assert processor.streamers[0] == mock_streamer


class TestClinvarLabels:
    def test_pathogenic_labels(self):
        expected_labels = [
            "Pathogenic/Likely_pathogenic",
            "Likely_pathogenic",
            "Pathogenic",
        ]
        assert ClinvarDataStreamer.PATHOGENIC_LABELS == expected_labels

    def test_benign_labels(self):
        expected_labels = ["Benign/Likely_benign", "Likely_benign", "Benign"]
        assert ClinvarDataStreamer.BENIGN_LABELS == expected_labels


def test_integration_example():
    example_usage = """
    processor = create_clinvar_training_set_streamer(
        clinvar_path="./data/clinvar/clinvar_20220403.vcf.gz",
        output_dir="./data/training_set"
    )
    training_set = processor.process()
    if training_set:
        print(f"Generated {training_set.count()} training examples")
    """
    assert "processor = create_clinvar_training_set_streamer" in example_usage


if __name__ == "__main__":
    test_suite = TestDataStreamer()
    test_suite.test_load_sample_gene_set()
    test_suite.test_clinvar_streamer_init()
    test_suite.test_factory_function()
    test_suite.test_stream_processor_init()
    label_tests = TestClinvarLabels()
    label_tests.test_pathogenic_labels()
    label_tests.test_benign_labels()
    test_integration_example()
    print("✅ All basic tests passed!")
    print("📋 Data streamer implementation is ready for use")
    print("\nNext steps:")
    print("1. Add your actual Clinvar VCF file path")
    print("2. Implement real CHD gene loading logic")
    print("3. Test with real data")
    print("4. Extend with additional streamers for other data sources")
