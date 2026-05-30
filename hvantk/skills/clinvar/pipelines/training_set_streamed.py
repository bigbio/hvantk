# Refactored Generate Training Set Script using Data Streamers
# eam
# 07.04.22
# Refactored to use streaming architecture

"""
Generate training set from Clinvar using data streamers
"""

import logging
import os
from hvantk.algorithms.hgc.constants import VCF_EXTENSION
from hvantk.core.utils.clinvar_streamer import create_clinvar_training_set_streamer
from hvantk.core.utils.gene_sets import load_sample_chd_gene_set

logger = logging.getLogger(__name__)


def main():
    """Main function to generate training set using data streamers"""

    # Configuration
    output_dir = "./data/training_set"
    clinvar_path = os.environ.get(
        "CLINVAR_VCF", f"./data/clinvar/clinvar_20220403{VCF_EXTENSION}"
    )

    logger.info(f"Starting Clinvar training set generation")
    logger.info(f"Clinvar path: {clinvar_path}")
    logger.info(f"Output directory: {output_dir}")

    gene_set = load_sample_chd_gene_set()  # default sandbox set
    logger.info(f"Loaded {len(gene_set)} sample genes")

    # Create the streaming processor
    processor = create_clinvar_training_set_streamer(
        clinvar_path=clinvar_path, output_dir=output_dir, gene_set=gene_set
    )

    # Process the data
    try:
        training_set = processor.process()
        if training_set:
            logger.info(
                f"Successfully generated training set with {training_set.count()} variants"
            )

            # Show some statistics
            tp_count = training_set.filter(training_set.rf_label == "TP").count()
            tn_count = training_set.filter(training_set.rf_label == "TN").count()

            logger.info(f"Training set statistics:")
            logger.info(f"  True Positives (TP): {tp_count}")
            logger.info(f"  True Negatives (TN): {tn_count}")
            logger.info(f"  Total: {tp_count + tn_count}")

        else:
            logger.warning("No training set generated")

    except Exception as e:
        logger.error(f"Error generating training set: {e}")
        raise


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    main()
