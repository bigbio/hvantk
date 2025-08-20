# Refactored Generate Training Set Script using Data Streamers
# eam
# 07.04.22
# Refactored to use streaming architecture

"""
Generate training set from Clinvar using data streamers
"""

import logging
from hvantk.utils.clinvar_streamer import create_clinvar_training_set_streamer, load_chd_gene_set
from hvantk.core.config import RAW_DATA_PATHS

logger = logging.getLogger(__name__)

def main():
    """Main function to generate training set using data streamers"""

    # Configuration
    output_dir = "./data/training_set"
    clinvar_path = RAW_DATA_PATHS.get("clinvar_path", "./data/clinvar/clinvar_20220403.vcf.gz")

    logger.info(f"Starting Clinvar training set generation")
    logger.info(f"Clinvar path: {clinvar_path}")
    logger.info(f"Output directory: {output_dir}")

    # Load CHD gene set
    logger.info("Loading CHD-associated genes")
    chd_gene_set = load_chd_gene_set()
    logger.info(f"Loaded {len(chd_gene_set)} CHD genes")

    # Create the streaming processor
    processor = create_clinvar_training_set_streamer(
        clinvar_path=clinvar_path,
        output_dir=output_dir,
        chd_genes=chd_gene_set
    )

    # Process the data
    try:
        training_set = processor.process()
        if training_set:
            logger.info(f"Successfully generated training set with {training_set.count()} variants")

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
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    main()
