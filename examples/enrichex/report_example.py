#!/usr/bin/env python3
"""
EnrichEx Report Example

Generates a self-contained HTML report from pre-computed EnrichEx results.
"""

import argparse
import logging
from pathlib import Path

from hvantk.algorithms.enrichex.report import generate_report

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Generate an EnrichEx HTML report.")
    parser.add_argument(
        "--output",
        default="reports/synthetic_enrichex_report.html",
        help="Output HTML path (default: reports/synthetic_enrichex_report.html)",
    )
    args = parser.parse_args()

    results_dir = Path(__file__).parent / "results"
    overlap_path = results_dir / "overlap_results.tsv"
    gene_sets_path = Path(__file__).parent / "synthetic_gene_sets.json"

    if not overlap_path.exists():
        logger.error(
            f"Results file not found: {overlap_path}\n"
            "Run 'hvantk enrichex overlap ...' or overlap_enrichment_example.py first."
        )
        return

    generate_report(
        output_path=args.output,
        overlap_results=str(overlap_path),
        gene_sets_path=str(gene_sets_path) if gene_sets_path.exists() else None,
        title="Synthetic EnrichEx Report",
        description="Example report generated from synthetic overlap enrichment results.",
    )
    logger.info(f"Report saved to {args.output}")


if __name__ == "__main__":
    main()
