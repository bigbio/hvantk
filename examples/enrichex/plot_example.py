#!/usr/bin/env python3
"""
EnrichEx Plot Example

Demonstrates how to generate publication-quality plots from pre-computed
EnrichEx results using the plotting API.
"""

import logging
from pathlib import Path

import pandas as pd

from hvantk.enrichex.plot import plot_enrichment_dotplot, plot_enrichment_barplot

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def main():
    results_dir = Path(__file__).parent / "results"
    output_dir = Path("plots/synthetic_enrichex")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load overlap results
    overlap_path = results_dir / "overlap_results.tsv"
    if not overlap_path.exists():
        logger.error(
            f"Results file not found: {overlap_path}\n"
            "Run 'hvantk enrichex overlap ...' or overlap_enrichment_example.py first."
        )
        return

    df = pd.read_csv(str(overlap_path), sep="\t")
    logger.info(f"Loaded {len(df)} enrichment results from {overlap_path}")

    # Generate dot plot
    dotplot_path = str(output_dir / "enrichex_dotplot.png")
    plot_enrichment_dotplot(df, output_path=dotplot_path, top_n=20)
    logger.info(f"Dot plot saved to {dotplot_path}")

    # Generate bar plot
    barplot_path = str(output_dir / "enrichex_barplot.png")
    plot_enrichment_barplot(df, output_path=barplot_path, top_n=20)
    logger.info(f"Bar plot saved to {barplot_path}")

    logger.info("Done! All plots saved to %s", output_dir)


if __name__ == "__main__":
    main()
