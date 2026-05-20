"""Hail Table builder for the GenCC (Gene Curation Coalition) submissions resource.

This module owns ``create_gencc_submissions_tb``, the canonical builder
that turns the GenCC submissions TSV into a Hail Table. It was migrated
out of :mod:`hvantk.core.builders.table` so that everything
GenCC-specific (builder, downloader, dataset class, tests, fixtures,
SKILL) lives under the plugin folder at :mod:`hvantk.skills.gencc`.

The shared helper ``_create_table_base`` and the ``get_row_fields``
utility intentionally stay in their existing modules because they are
reused by other builders.
"""

from __future__ import annotations

import logging
from typing import List, Optional

import hail as hl

from hvantk.core.constants import (
    GENCC_CLASSIFICATION_LEVELS,
    GENCC_SUBMISSION_FIELDS,
)
from hvantk.core.builders.table import _create_table_base
from hvantk.core.utils.table_utils import get_row_fields

logger = logging.getLogger(__name__)


def create_gencc_submissions_tb(
    input_path: str,
    output_path: str,
    key_by: str = "gene_disease_submitter",
    min_classification: Optional[str] = None,
    fields: Optional[List[str]] = None,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> "hl.Table":
    """
    Create a Hail Table from a GenCC submissions TSV file.

    GenCC aggregates gene-disease validity assertions from 12+ submitting
    organizations (ClinGen, PanelApp, G2P, Orphanet, etc.).

    Parameters
    ----------
    input_path : str
        Path to the GenCC submissions TSV file.
    output_path : str
        Path to write the output Hail Table.
    key_by : str, optional
        Keying strategy:
        - "gene_disease_submitter" (default): Key by (hgnc_id, mondo_id, submitter)
        - "gene_disease": Aggregate across submitters, key by (hgnc_id, mondo_id)
        - "gene": Aggregate all diseases per gene, key by hgnc_id
    min_classification : str, optional
        Filter to classifications at or above this level.
    fields : list of str, optional
        List of fields to select from the table (default: None, keeps all).
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).

    Returns
    -------
    hl.Table
        Hail Table with GenCC gene-disease-submitter validity annotations.
    """
    valid_keys = ("gene_disease_submitter", "gene_disease", "gene")
    if key_by not in valid_keys:
        raise ValueError(f"key_by must be one of {valid_keys}, got: {key_by}")

    if (
        min_classification is not None
        and min_classification not in GENCC_CLASSIFICATION_LEVELS
    ):
        raise ValueError(
            f"min_classification must be one of {GENCC_CLASSIFICATION_LEVELS}, "
            f"got: {min_classification}"
        )

    def import_func():
        return hl.import_table(
            paths=input_path,
            delimiter="\t",
            impute=False,
            min_partitions=10,
        )

    def transform(ht: hl.Table) -> hl.Table:
        # Rename fields to standardized names
        logger.info("Renaming GenCC fields to standardized names")
        rename_map = {
            k: v for k, v in GENCC_SUBMISSION_FIELDS.items() if k in get_row_fields(ht)
        }
        ht = ht.rename(rename_map)

        # Clean HGNC ID (strip "HGNC:" prefix)
        row_fields = get_row_fields(ht)
        if "hgnc_id" in row_fields:
            ht = ht.annotate(
                hgnc_id=hl.if_else(
                    ht.hgnc_id.startswith("HGNC:"),
                    ht.hgnc_id.replace("HGNC:", ""),
                    ht.hgnc_id,
                )
            )

        # Clean MONDO ID (strip "MONDO:" prefix if present)
        if "mondo_id" in row_fields:
            ht = ht.annotate(
                mondo_id=hl.if_else(
                    ht.mondo_id.startswith("MONDO:"),
                    ht.mondo_id.replace("MONDO:", ""),
                    ht.mondo_id,
                )
            )

        # Add classification level as numeric for filtering/sorting
        classification_order = {
            level: i for i, level in enumerate(GENCC_CLASSIFICATION_LEVELS)
        }
        ht = ht.annotate(
            classification_level=hl.literal(classification_order).get(
                ht.classification, hl.len(GENCC_CLASSIFICATION_LEVELS)
            )
        )

        # Apply min_classification filter if specified
        if min_classification is not None:
            min_level = classification_order[min_classification]
            logger.info(
                f"Filtering to classifications >= {min_classification} "
                f"(level {min_level})"
            )
            ht = ht.filter(ht.classification_level <= min_level)

        # Apply keying strategy
        if key_by == "gene_disease_submitter":
            logger.info("Keying by (hgnc_id, mondo_id, submitter)")
            ht = ht.key_by("hgnc_id", "mondo_id", "submitter")
        elif key_by == "gene_disease":
            logger.info("Aggregating by gene-disease (hgnc_id, mondo_id)")
            ht = (
                ht.group_by("hgnc_id", "mondo_id", "gene_symbol", "disease_label")
                .aggregate(
                    submitters=hl.agg.collect_as_set(ht.submitter),
                    classifications=hl.agg.collect_as_set(ht.classification),
                    modes_of_inheritance=hl.agg.collect_as_set(ht.mode_of_inheritance),
                    max_classification_level=hl.agg.min(ht.classification_level),
                )
                .key_by("hgnc_id", "mondo_id")
            )
            ht = ht.annotate(n_submitters=hl.len(ht.submitters))
            classification_labels = hl.literal(GENCC_CLASSIFICATION_LEVELS)
            safe_index = hl.min(
                ht.max_classification_level, hl.len(classification_labels) - 1
            )
            ht = ht.annotate(
                max_classification_label=classification_labels[safe_index],
                classification=classification_labels[safe_index],
                classification_level=ht.max_classification_level,
            )
        else:  # key_by == "gene"
            logger.info("Aggregating by gene (hgnc_id)")
            ht = (
                ht.group_by("hgnc_id", "gene_symbol")
                .aggregate(
                    disease_labels=hl.agg.collect_as_set(ht.disease_label),
                    disease_mondo_pairs=hl.agg.collect_as_set(
                        hl.struct(disease_label=ht.disease_label, mondo_id=ht.mondo_id)
                    ),
                    mondo_ids=hl.agg.collect_as_set(ht.mondo_id),
                    classifications=hl.agg.collect_as_set(ht.classification),
                    modes_of_inheritance=hl.agg.collect_as_set(ht.mode_of_inheritance),
                    submitters=hl.agg.collect_as_set(ht.submitter),
                    max_classification_level=hl.agg.min(ht.classification_level),
                )
                .key_by("hgnc_id")
            )
            ht = ht.annotate(
                n_diseases=hl.len(ht.mondo_ids),
                n_submitters=hl.len(ht.submitters),
            )
            classification_labels = hl.literal(GENCC_CLASSIFICATION_LEVELS)
            safe_index = hl.min(
                ht.max_classification_level, hl.len(classification_labels) - 1
            )
            ht = ht.annotate(max_classification_label=classification_labels[safe_index])

        return ht

    gencc_tb = _create_table_base(
        source_name="GenCC Submissions",
        input_path=input_path,
        output_path=output_path,
        import_func=import_func,
        transform_func=transform,
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )

    return gencc_tb


def build_gencc_submissions(
    parsed_input,
    ctx,
    *,
    key_by: str = "gene_disease_submitter",
    min_classification=None,
    fields=None,
):
    """Phase B builder — returns an AnnotationTable.

    See create_gencc_submissions_tb for the transform semantics; this function
    runs the same import + transform but returns the lazy Hail Table wrapped
    in an AnnotationTable with Provenance.
    """
    from hvantk.core.models import AnnotationTable

    valid_keys = ("gene_disease_submitter", "gene_disease", "gene")
    if key_by not in valid_keys:
        raise ValueError(f"key_by must be one of {valid_keys}, got: {key_by}")
    if (
        min_classification is not None
        and min_classification not in GENCC_CLASSIFICATION_LEVELS
    ):
        raise ValueError(
            f"min_classification must be one of {GENCC_CLASSIFICATION_LEVELS}, "
            f"got: {min_classification}"
        )

    ht = hl.import_table(
        paths=str(parsed_input),
        delimiter="\t",
        impute=False,
        min_partitions=10,
    )

    rename_map = {
        k: v for k, v in GENCC_SUBMISSION_FIELDS.items() if k in get_row_fields(ht)
    }
    ht = ht.rename(rename_map)

    row_fields = get_row_fields(ht)
    if "hgnc_id" in row_fields:
        ht = ht.annotate(
            hgnc_id=hl.if_else(
                ht.hgnc_id.startswith("HGNC:"),
                ht.hgnc_id.replace("HGNC:", ""),
                ht.hgnc_id,
            )
        )
    if "mondo_id" in row_fields:
        ht = ht.annotate(
            mondo_id=hl.if_else(
                ht.mondo_id.startswith("MONDO:"),
                ht.mondo_id.replace("MONDO:", ""),
                ht.mondo_id,
            )
        )

    classification_order = {
        level: i for i, level in enumerate(GENCC_CLASSIFICATION_LEVELS)
    }
    ht = ht.annotate(
        classification_level=hl.literal(classification_order).get(
            ht.classification, hl.len(GENCC_CLASSIFICATION_LEVELS)
        )
    )

    if min_classification is not None:
        min_level = classification_order[min_classification]
        ht = ht.filter(ht.classification_level <= min_level)

    # Apply default keying (gene_disease_submitter); leave aggregation modes
    # to the legacy function for now.
    if key_by != "gene_disease_submitter":
        raise NotImplementedError(
            f"Phase B builder currently supports key_by='gene_disease_submitter' only; "
            f"got {key_by}. Use the legacy create_gencc_submissions_tb function for now."
        )
    ht = ht.key_by("hgnc_id", "mondo_id", "submitter")

    if fields is not None:
        ht = ht.select(*fields)

    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="gencc-submissions-v1")
    )
