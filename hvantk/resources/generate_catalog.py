"""
Master catalog for the reorganized hvantk registry system.
Provides unified access to all omics datasets.
"""

import json
import yaml
from pathlib import Path
from typing import Dict, List, Any
from datetime import datetime


def generate_master_catalog():
    """Generate a master catalog.yaml file with all datasets from the new registry structure."""

    registry_root = Path(__file__).parent / "registry"
    omics_types = ["transcriptomics", "proteomics", "genomics", "epigenomics"]

    catalog = {
        "version": "2.0",
        "generated": datetime.now().isoformat(),
        "description": "Unified catalog of all omics datasets in hvantk registry",
        "registry_structure": {
            "schemas": "hvantk/resources/schemas/",
            "transcriptomics": "hvantk/resources/registry/transcriptomics/",
            "proteomics": "hvantk/resources/registry/proteomics/",
            "genomics": "hvantk/resources/registry/genomics/",
            "epigenomics": "hvantk/resources/registry/epigenomics/",
        },
        "datasets": {},
    }

    total_datasets = 0

    for omics_type in omics_types:
        datasets_file = registry_root / omics_type / "datasets.json"

        if datasets_file.exists():
            with open(datasets_file, "r") as f:
                datasets = json.load(f)

            catalog["datasets"][omics_type] = {
                "count": len(datasets),
                "schemas": f"schemas/{omics_type}_schema.json",
                "data_file": f"registry/{omics_type}/datasets.json",
                "summary": {
                    "data_sources": list(
                        set(d.get("data_source", "Unknown") for d in datasets)
                    ),
                    "organisms": list(
                        set(d.get("organism", "Unknown") for d in datasets)
                    ),
                    "sample_range": {
                        "min": min(
                            (d.get("sample_count", 0) for d in datasets), default=0
                        ),
                        "max": max(
                            (d.get("sample_count", 0) for d in datasets), default=0
                        ),
                    },
                },
            }

            total_datasets += len(datasets)
        else:
            catalog["datasets"][omics_type] = {
                "count": 0,
                "schemas": f"schemas/{omics_type}_schema.json",
                "data_file": f"registry/{omics_type}/datasets.json",
                "summary": {
                    "data_sources": [],
                    "organisms": [],
                    "sample_range": {"min": 0, "max": 0},
                },
            }

    catalog["summary"] = {
        "total_datasets": total_datasets,
        "omics_types": len(
            [t for t in omics_types if catalog["datasets"][t]["count"] > 0]
        ),
    }

    # Save catalog
    catalog_file = Path(__file__).parent / "catalog.yaml"
    with open(catalog_file, "w") as f:
        yaml.dump(catalog, f, default_flow_style=False, sort_keys=False)

    print(f"✓ Generated master catalog with {total_datasets} total datasets")
    return catalog


if __name__ == "__main__":
    generate_master_catalog()
