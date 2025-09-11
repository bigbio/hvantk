#!/usr/bin/env python3
"""Create status badge for dataset validation."""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path.cwd()))

try:
    from hvantk.register import RegistryManager

    manager = RegistryManager("./dataset_validation_results/validation_registry.json")
    stats = manager.get_registry_statistics() if manager else {}

    rate = stats.get("success_rate", 0)
    successful = stats.get("successful", 0)
    total = stats.get("total_datasets", 0)

    print(f"Success Rate: {rate}% ({successful}/{total})")

    # Ensure output directory exists
    Path("./docs/registry").mkdir(parents=True, exist_ok=True)

    # Write status to file
    with open("./docs/registry/status.txt", "w") as f:
        f.write(str(rate))

    print("✅ Status badge created successfully")

except Exception as e:
    print(f"❌ Error creating status badge: {e}")
    # Create minimal status file
    Path("./docs/registry").mkdir(parents=True, exist_ok=True)
    with open("./docs/registry/status.txt", "w") as f:
        f.write("0")
