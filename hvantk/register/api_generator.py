#!/usr/bin/env python3
"""
API Endpoint Generator for HVANTK Dataset Validation Registry

This module generates RESTful API endpoints for programmatic access to
dataset validation registry data.
"""
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any
from collections import Counter

from .config import RegistryConfig


class APIEndpointGenerator:
    """Generator for API endpoints from validation registry data."""

    def __init__(self, config: RegistryConfig = None):
        """Initialize the API generator with configuration."""
        self.config = config or RegistryConfig()

    def generate_status_api(self, registry: Dict[str, Any]) -> Dict[str, Any]:
        """Generate /api/status.json endpoint."""
        total_datasets = len(registry)
        status_counts = Counter(data['status'] for data in registry.values())

        successful = sum(count for status, count in status_counts.items()
                        if status in ['tier3_passed', 'tier2_passed'])

        return {
            "status": "ok",
            "timestamp": datetime.now().isoformat(),
            "summary": {
                "total_datasets": total_datasets,
                "successful_validations": successful,
                "failed_validations": total_datasets - successful,
                "success_rate": round((successful / total_datasets * 100) if total_datasets > 0 else 0, 1)
            },
            "status_breakdown": dict(status_counts),
            "last_update": max((data['timestamp'] for data in registry.values()), default=None)
        }

    def generate_datasets_api(self, registry: Dict[str, Any]) -> Dict[str, Any]:
        """Generate /api/datasets.json endpoint."""
        return {
            "datasets": registry,
            "metadata": {
                "total_count": len(registry),
                "generated_at": datetime.now().isoformat(),
                "schema_version": "1.0"
            }
        }

    def generate_stats_api(self, registry: Dict[str, Any]) -> Dict[str, Any]:
        """Generate /api/stats.json endpoint."""
        # Calculate statistics
        status_counts = Counter(data['status'] for data in registry.values())
        type_counts = Counter(data['dataset_type'] for data in registry.values())

        # Time-based analysis
        recent_failures = []
        for dataset_id, data in registry.items():
            if data['status'] in ['tier1_failed', 'tier2_failed', 'tier3_failed', 'download_failed']:
                recent_failures.append({
                    "dataset_id": dataset_id,
                    "status": data['status'],
                    "timestamp": data['timestamp'],
                    "error_message": data.get('error_message', '')[:200]  # Truncate long errors
                })

        # Sort by timestamp (most recent first)
        recent_failures.sort(key=lambda x: x['timestamp'], reverse=True)

        return {
            "statistics": {
                "by_status": dict(status_counts),
                "by_type": dict(type_counts),
                "success_rates": {
                    "overall": round((sum(count for status, count in status_counts.items()
                                        if status in ['tier3_passed', 'tier2_passed']) / len(registry) * 100)
                                   if len(registry) > 0 else 0, 1),
                    "ucsc": self.calculate_success_rate(registry, 'ucsc'),
                    "expression_atlas": self.calculate_success_rate(registry, 'expression_atlas')
                }
            },
            "recent_failures": recent_failures[:10],  # Last 10 failures
            "generated_at": datetime.now().isoformat()
        }

    def calculate_success_rate(self, registry: Dict[str, Any], dataset_type: str) -> float:
        """Calculate success rate for a specific dataset type."""
        type_datasets = [data for data in registry.values() if data['dataset_type'] == dataset_type]
        if not type_datasets:
            return 0.0

        successful = sum(1 for data in type_datasets
                        if data['status'] in ['tier3_passed', 'tier2_passed'])
        return round((successful / len(type_datasets) * 100), 1)

    def generate_api_documentation(self) -> Dict[str, Any]:
        """Generate API documentation."""
        return {
            "api_version": "1.0",
            "endpoints": {
                "/api/status.json": {
                    "description": "Overall status and summary statistics",
                    "methods": ["GET"],
                    "response_format": "JSON"
                },
                "/api/datasets.json": {
                    "description": "Complete dataset validation data",
                    "methods": ["GET"],
                    "response_format": "JSON"
                },
                "/api/stats.json": {
                    "description": "Detailed statistics and recent failures",
                    "methods": ["GET"],
                    "response_format": "JSON"
                }
            },
            "usage_examples": {
                "curl": "curl https://username.github.io/repo/api/status.json",
                "python": "import requests; r = requests.get('https://username.github.io/repo/api/status.json')",
                "javascript": "fetch('https://username.github.io/repo/api/status.json').then(r => r.json())"
            },
            "last_updated": datetime.now().isoformat()
        }

    def generate_api_endpoints(self, registry_file: str, output_dir: str) -> Path:
        """Generate all API endpoints from validation registry."""
        # Load registry data
        with open(registry_file, 'r') as f:
            registry = json.load(f)

        # Create API output directory
        api_dir = Path(output_dir) / 'api'
        api_dir.mkdir(parents=True, exist_ok=True)

        # Generate API endpoints
        endpoints = {
            'status.json': self.generate_status_api(registry),
            'datasets.json': self.generate_datasets_api(registry),
            'stats.json': self.generate_stats_api(registry)
        }

        # Write API files
        for filename, data in endpoints.items():
            with open(api_dir / filename, 'w') as f:
                json.dump(data, f, indent=2)
            print(f"✅ Generated {filename}")

        # Generate API documentation
        docs = self.generate_api_documentation()
        with open(api_dir / 'index.json', 'w') as f:
            json.dump(docs, f, indent=2)

        print(f"🌐 API documentation generated at {api_dir / 'index.json'}")
        print(f"📊 Total API endpoints: {len(endpoints)}")

        return api_dir
