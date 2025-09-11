#!/usr/bin/env python3
"""Generate validation summary for dataset validation."""

import sys
import os
from pathlib import Path
from datetime import datetime

# Add project root to path
sys.path.insert(0, str(Path.cwd()))

try:
    from hvantk.register import RegistryManager
    manager = RegistryManager("./dataset_validation_results/validation_registry.json")
    stats = manager.get_registry_statistics()
    failed_datasets = manager.get_failed_datasets(5)
except Exception as e:
    print(f"Warning: Could not load registry data: {e}")
    stats = {"total_datasets": 0, "successful": 0, "failed": 0, "success_rate": 0.0}
    failed_datasets = []

run_id = os.environ.get("GITHUB_RUN_NUMBER", "unknown")
timestamp = datetime.utcnow().isoformat() + "Z"

lines = []
lines.append("## Dataset Validation Summary 📊")
lines.append("")
lines.append(f"**Run ID:** {run_id}")
lines.append(f"**Timestamp:** {timestamp}")
lines.append("")
lines.append("### Overall Results")
lines.append(f"- 📈 **Total Datasets:** {stats.get('total_datasets', 0)}")
lines.append(f"- ✅ **Successful:** {stats.get('successful', 0)} ({stats.get('success_rate', 0)}%)")
lines.append(f"- ❌ **Failed:** {stats.get('failed', 0)}")

by_status = stats.get("by_status") or {}
if by_status:
    lines.append("")
    lines.append("### Status Breakdown")
    for status, count in sorted(by_status.items()):
        emoji = "✅" if "passed" in status else ("❌" if "failed" in status else "⚪")
        status_title = status.replace("_", " ").title()
        lines.append(f"- {emoji} **{status_title}:** {count}")

if failed_datasets:
    lines.append("")
    lines.append(f"### Recent Failures ({len(failed_datasets)})")
    for failure in failed_datasets:
        msg = (failure.get("error_message") or "")
        err = (msg[:100] + "...") if len(msg) > 100 else msg
        lines.append(f"- **{failure.get('dataset_id','unknown')}** ({failure.get('status','unknown')}): {err}")

with open("validation_summary.md", "w") as f:
    f.write("\n".join(lines))

print("✅ Validation summary generated")
