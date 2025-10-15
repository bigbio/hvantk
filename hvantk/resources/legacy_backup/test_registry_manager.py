"""
Test the new registry management framework.
"""
import sys
from pathlib import Path

# Add the parent directory to path
sys.path.append(str(Path(__file__).parent))
from registry_manager import RegistryManager

def test_registry_manager():
    """Test the new registry management system."""
    print("Testing Registry Manager...")

    # Initialize the registry manager
    registry = RegistryManager()

    # Get registry statistics
    stats = registry.get_registry_stats()
    print(f"\n=== Registry Statistics ===")
    print(f"Total datasets: {stats['total_datasets']}")
    print(f"By omics type: {stats['by_omics_type']}")
    print(f"By organism: {dict(list(stats['by_organism'].items())[:5])}")  # Show first 5
    print(f"By data source: {stats['by_data_source']}")

    # Test search functionality
    print(f"\n=== Search Tests ===")

    # Search for GTEx datasets
    gtex_results = registry.search_datasets("GTEx")
    print(f"GTEx search results: {len(gtex_results)} datasets")

    # Search for brain datasets
    brain_results = registry.search_datasets("brain")
    print(f"Brain search results: {len(brain_results)} datasets")

    # Search by organism
    mouse_results = registry.search_datasets(organism="Mus musculus")
    print(f"Mouse datasets: {len(mouse_results)}")

    # Test getting specific dataset
    gtex_dataset = registry.get_dataset_by_accession("E-GTEX-8")
    if gtex_dataset:
        print(f"\nFound GTEx dataset: {gtex_dataset['title']}")
        print(f"Omics type: {gtex_dataset['_omics_type']}")

    print("\n✓ Registry Manager tests completed successfully!")

if __name__ == "__main__":
    test_registry_manager()
