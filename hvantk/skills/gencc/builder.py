"""Hail Table builder for the GenCC (Gene Curation Coalition) submissions resource.

Owns the Phase B ``build_gencc_submissions`` builder. GenCC aggregates
gene-disease validity assertions from 12+ submitting organizations (ClinGen,
PanelApp, G2P, Orphanet, ...). Keys by ``(hgnc_id, mondo_id, submitter)``.
"""

from __future__ import annotations

import logging

import hail as hl

from hvantk.core.constants import (
    GENCC_CLASSIFICATION_LEVELS,
    GENCC_SUBMISSION_FIELDS,
)
from hvantk.core.utils.table_utils import get_row_fields

logger = logging.getLogger(__name__)


def build_gencc_submissions(
    parsed_input,
    ctx,
    *,
    min_classification=None,
    fields=None,
):
    """Phase B builder — returns an AnnotationTable keyed by
    ``(hgnc_id, mondo_id, submitter)``.

    Optional ``min_classification`` filters to assertions at or above the
    given confidence level (one of ``GENCC_CLASSIFICATION_LEVELS``).
    """
    from hvantk.core.models import AnnotationTable

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

    ht = ht.key_by("hgnc_id", "mondo_id", "submitter")

    if fields is not None:
        ht = ht.select(*fields)

    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="gencc-submissions-v1")
    )
