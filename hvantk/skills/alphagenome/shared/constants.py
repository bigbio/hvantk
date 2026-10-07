"""Plugin-local constants for the alphagenome skill."""

#: Exact ``variant_scorer`` string written by the AlphaGenome SDK's
#: ``variant_scorers.tidy_scores()`` -> the summary field the builder emits for it.
#:
#: These are the 19 scorers of the SDK's ``RECOMMENDED_VARIANT_SCORERS``, copied
#: verbatim from a real ``tidy_scores()`` shard (16.6M rows, 411 ClinVar variants,
#: inspected 2026-10-06). Each scorer reports exactly one ``output_type``, so the
#: field name carries it; the builder reads ``output_type`` from its own column and
#: never parses it from these strings, because three of them (``PolyadenylationScorer()``,
#: ``ContactMapScorer()``, ``SpliceJunctionScorer()``) do not name it.
#:
#: A scorer missing from this mapping makes the builder raise rather than guess: a new
#: SDK release can add scorers or change how they print (see SKILL.md s 8).
SCORER_FIELDS: dict[str, str] = {
    # ATAC
    "CenterMaskScorer(requested_output=ATAC, width=501, aggregation_type=ACTIVE_SUM)": "atac_active_sum",
    "CenterMaskScorer(requested_output=ATAC, width=501, aggregation_type=DIFF_LOG2_SUM)": "atac_diff_log2_sum",
    # CAGE
    "CenterMaskScorer(requested_output=CAGE, width=501, aggregation_type=ACTIVE_SUM)": "cage_active_sum",
    "CenterMaskScorer(requested_output=CAGE, width=501, aggregation_type=DIFF_LOG2_SUM)": "cage_diff_log2_sum",
    # CHIP_HISTONE
    "CenterMaskScorer(requested_output=CHIP_HISTONE, width=2001, aggregation_type=ACTIVE_SUM)": "chip_histone_active_sum",
    "CenterMaskScorer(requested_output=CHIP_HISTONE, width=2001, aggregation_type=DIFF_LOG2_SUM)": (
        "chip_histone_diff_log2_sum"
    ),
    # CHIP_TF
    "CenterMaskScorer(requested_output=CHIP_TF, width=501, aggregation_type=ACTIVE_SUM)": "chip_tf_active_sum",
    "CenterMaskScorer(requested_output=CHIP_TF, width=501, aggregation_type=DIFF_LOG2_SUM)": "chip_tf_diff_log2_sum",
    # CONTACT_MAPS
    "ContactMapScorer()": "contact_maps",
    # DNASE
    "CenterMaskScorer(requested_output=DNASE, width=501, aggregation_type=ACTIVE_SUM)": "dnase_active_sum",
    "CenterMaskScorer(requested_output=DNASE, width=501, aggregation_type=DIFF_LOG2_SUM)": "dnase_diff_log2_sum",
    # PROCAP
    "CenterMaskScorer(requested_output=PROCAP, width=501, aggregation_type=ACTIVE_SUM)": "procap_active_sum",
    "CenterMaskScorer(requested_output=PROCAP, width=501, aggregation_type=DIFF_LOG2_SUM)": "procap_diff_log2_sum",
    # RNA_SEQ
    "GeneMaskActiveScorer(requested_output=RNA_SEQ)": "rna_seq_active",
    "GeneMaskLFCScorer(requested_output=RNA_SEQ)": "rna_seq_lfc",
    "PolyadenylationScorer()": "rna_seq_polyadenylation",
    # SPLICE_JUNCTIONS
    "SpliceJunctionScorer()": "splice_junctions",
    # SPLICE_SITE_USAGE
    "GeneMaskSplicingScorer(requested_output=SPLICE_SITE_USAGE, width=None)": "splice_site_usage",
    # SPLICE_SITES
    "GeneMaskSplicingScorer(requested_output=SPLICE_SITES, width=None)": "splice_sites",
}

#: The ``output_type`` values those scorers report; the ``output_types`` builder
#: parameter is checked against this set.
OUTPUT_TYPES: frozenset[str] = frozenset(
    {
        "ATAC",
        "CAGE",
        "CHIP_HISTONE",
        "CHIP_TF",
        "CONTACT_MAPS",
        "DNASE",
        "PROCAP",
        "RNA_SEQ",
        "SPLICE_JUNCTIONS",
        "SPLICE_SITE_USAGE",
        "SPLICE_SITES",
    }
)
