#!/usr/bin/env python3
"""
CHD cell-type burden analysis — ready for external host.

Adapts BurdenPipeline to the CHD cohort MatrixTable schema:
  - Consequence field: 'Consequence' (capital C)
  - AF field: 'gnomAD_AF' (camelCase)
  - Prediction score field: 'vep.CADD_PHRED' (nested, capital)
  - Phenotype: 'phe.is_case' in MT columns (extracted automatically)
  - No GQ/DP in entries, no filters field → disable QC filters

Usage:
    python chd_burden_run.py \
        --cohort-mt /path/to/cohort.mt \
        --gene-sets heart:/path/to/heart_cell_types.json \
        --gene-sets brain:/path/to/brain_cell_types.json \
        --output-dir /path/to/results \
        [--dry-run]
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import Dict

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(name)s [%(levelname)s] %(message)s",
    stream=sys.stdout,
)
logger = logging.getLogger("chd_burden")


def parse_args():
    p = argparse.ArgumentParser(description="CHD cell-type burden analysis")
    p.add_argument("--cohort-mt", required=True, help="Path to cohort MatrixTable")
    p.add_argument(
        "--gene-sets",
        action="append",
        required=True,
        help="name:/path/to/collection.json (repeatable)",
    )
    p.add_argument(
        "--output-dir", default="results/chd_burden", help="Output directory"
    )
    p.add_argument(
        "--phenotype-field",
        default="phe.is_case",
        help="Phenotype column field (dot notation for nested structs)",
    )
    p.add_argument(
        "--covariates",
        default=None,
        help="Comma-separated covariate fields from MT cols (dot notation ok)",
    )
    p.add_argument(
        "--max-af", type=float, default=0.001, help="Max gnomAD AF (default: 0.001)"
    )
    p.add_argument(
        "--min-cadd", type=float, default=None, help="Min CADD PHRED (optional)"
    )
    p.add_argument(
        "--min-carriers", type=int, default=5, help="Min carriers per gene set"
    )
    p.add_argument(
        "--competitive", action="store_true", help="Run permutation competitive test"
    )
    p.add_argument(
        "--n-permutations",
        type=int,
        default=10000,
        help="Permutations for competitive test",
    )
    p.add_argument("--dry-run", action="store_true", help="Show plan without running")
    p.add_argument(
        "--generate-report", action="store_true", help="Generate HTML report"
    )
    return p.parse_args()


def build_variant_classes(max_af: float, min_score: float = None) -> Dict:
    """Build variant class filters for the CHD schema.

    Maps to the CHD MT fields:
      - Consequence (capital C)
      - gnomAD_AF (camelCase)
      - vep.CADD_PHRED (nested)
    """
    from hvantk.enrichex.burden import VariantFilter

    # Common filter settings for this MT schema
    common = dict(
        af_field="gnomAD_AF",
        score_field="vep.CADD_PHRED",
        consequence_field="Consequence",
        max_af=max_af,
        pass_only=False,  # no 'filters' field in this MT
        min_gq=0,  # no GQ in entries
        min_dp=0,  # no DP in entries
    )

    classes = {
        "lof": VariantFilter(
            consequences=[
                "stop_gained",
                "frameshift_variant",
                "splice_donor_variant",
                "splice_acceptor_variant",
            ],
            min_score=None,  # LoF doesn't need CADD filtering
            **common,
        ),
        "missense_constrained": VariantFilter(
            consequences=["missense_variant"],
            min_score=min_score if min_score is not None else 25.0,
            **common,
        ),
        "synonymous": VariantFilter(
            consequences=["synonymous_variant"],
            min_score=None,  # negative control — no CADD filtering
            **common,
        ),
    }
    return classes


def main():
    args = parse_args()

    from hvantk.core.hail_context import init_hail

    init_hail(app_name="chd-burden")

    # Parse gene set collections
    gene_set_collections = {}
    for gs_arg in args.gene_sets:
        if ":" not in gs_arg:
            logger.error("Gene set arg must be name:/path — got: %s", gs_arg)
            sys.exit(1)
        name, path = gs_arg.split(":", 1)
        if not Path(path).exists():
            logger.error("Gene set file not found: %s", path)
            sys.exit(1)
        gene_set_collections[name] = path

    covariate_fields = args.covariates.split(",") if args.covariates else []

    # ── Build config ──────────────────────────────────────────────────
    # No need to manually extract phenotype — the pipeline handles
    # MT column fields natively when phenotype_ht_path is empty.
    from hvantk.enrichex.pipeline import BurdenConfig, BurdenPipeline

    config = BurdenConfig(
        cohort_mt_path=args.cohort_mt,
        # phenotype_ht_path omitted → extract from MT column fields
        phenotype_field=args.phenotype_field,  # e.g. "phe.is_case"
        phenotype_type="binary",
        covariate_fields=covariate_fields,
        gene_set_collections=gene_set_collections,
        variant_classes=build_variant_classes(
            max_af=args.max_af,
            min_score=args.min_score,
        ),
        gene_field="SYMBOL",
        min_carriers=args.min_carriers,
        correction_method="benjamini-hochberg",
        alpha=0.05,
        competitive=args.competitive,
        n_permutations=args.n_permutations,
        output_dir=args.output_dir,
        generate_report=args.generate_report,
    )

    pipeline = BurdenPipeline(config)
    pipeline.show_plan()

    if args.dry_run:
        logger.info("Dry run — exiting.")
        return

    # ── Run ────────────────────────────────────────────────────────────
    t0 = time.time()
    combined_df = pipeline.run()
    elapsed = time.time() - t0

    if combined_df.empty:
        logger.warning("No results produced.")
        return

    logger.info("")
    logger.info("=" * 70)
    logger.info("RESULTS SUMMARY")
    logger.info("=" * 70)
    logger.info("Total results: %d rows", len(combined_df))
    logger.info("Total time: %.1fs", elapsed)

    # Top hits
    sig = combined_df[combined_df.get("significant", False) == True]
    if not sig.empty:
        logger.info("\nSignificant results (%d):", len(sig))
        cols = ["gene_set_name", "variant_class", "collection", "p_value", "p_adjusted"]
        display_cols = [c for c in cols if c in sig.columns]
        print(sig.sort_values("p_adjusted")[display_cols].to_string(index=False))
    else:
        logger.info("No significant results at alpha=0.05.")

    # Synonymous control check
    syn = combined_df[combined_df["variant_class"] == "synonymous"]
    syn_sig = syn[syn.get("significant", False) == True] if not syn.empty else syn
    if not syn_sig.empty:
        logger.warning(
            "WARNING: %d gene sets significant for synonymous variants — "
            "possible confounding (gene length bias, population stratification).",
            len(syn_sig),
        )
        cols = ["gene_set_name", "collection", "p_adjusted"]
        display_cols = [c for c in cols if c in syn_sig.columns]
        print(
            syn_sig.sort_values("p_adjusted")[display_cols]
            .head(10)
            .to_string(index=False)
        )

    logger.info("\nOutputs in: %s", args.output_dir)


if __name__ == "__main__":
    main()
