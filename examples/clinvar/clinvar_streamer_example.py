# Example usage of hvantk's ClinVar variant streamer
# Demonstrates ClinVarVariantTableStreamer: loading a built ClinVar AnnotationTable,
# filtering by gene set and by pathogenicity, aggregating clinical-significance counts,
# and deriving a TP/TN training-label column.
#
# Build the Hail Table first (see README.md):
#
#     hvantk reprocess clinvar:variants --raw-dir data/clinvar --output clinvar.ht
#
# Usage:
#     python examples/clinvar/clinvar_streamer_example.py [path/to/clinvar.ht]

import logging
import sys

import hail as hl

from hvantk.skills.clinvar.streamers import ClinVarVariantTableStreamer

logger = logging.getLogger(__name__)

DEFAULT_TABLE_PATH = "clinvar.ht"


def main(table_path: str = DEFAULT_TABLE_PATH) -> None:
    hl.init(quiet=True)

    print("=== ClinVar Streamer Example ===")
    streamer = ClinVarVariantTableStreamer.from_path(table_path)
    ht = streamer.to_hail()
    print(f"Total variants: {ht.count()}")

    print("\n--- Filtering by gene set ---")
    gene_set = {"BRCA1", "BRCA2", "TP53"}
    panel_ht = streamer.filter_to_genes(gene_set)
    print(f"Variants in {sorted(gene_set)}: {panel_ht.count()}")

    print("\n--- Filtering by pathogenicity ---")
    pathogenic_ht = streamer.filter_by_pathogenicity(
        ["Pathogenic", "Likely_pathogenic"]
    )
    print(f"Pathogenic / likely pathogenic variants: {pathogenic_ht.count()}")

    print("\n--- Aggregating results ---")
    by_sig = panel_ht.aggregate(hl.agg.counter(panel_ht.info.CLNSIG))
    print(f"Clinical significance counts in the gene panel: {by_sig}")

    print("\n--- Deriving a training-label column ---")
    labeled_ht = streamer.label_training_set(
        gene_set={"BRCA1", "BRCA2"},
        label_column="rf_label",
    )
    labeled_ht = labeled_ht.filter(hl.is_defined(labeled_ht.rf_label))
    total_count = labeled_ht.count()
    tp_count = labeled_ht.filter(labeled_ht.rf_label == "TP").count()
    tn_count = total_count - tp_count
    print(f"Labeled training variants: {total_count}")
    if total_count:
        print(f"  True Positives: {tp_count} ({tp_count / total_count * 100:.1f}%)")
        print(f"  True Negatives: {tn_count} ({tn_count / total_count * 100:.1f}%)")

    hl.stop()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_TABLE_PATH)
