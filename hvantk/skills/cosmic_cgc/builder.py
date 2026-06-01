"""Hail Table builder for the COSMIC Cancer Gene Census (CGC) resource.

Phase B builder: imports the COSMIC CGC TSV, normalises tier classifications,
optionally filters by mutation context, and emits an AnnotationTable with
source-fingerprint provenance.
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Optional

import hail as hl

from hvantk.skills.cosmic_cgc.shared.constants import (
    COSMIC_CGC_FIELDS,
    COSMIC_CGC_CLASSIFICATION_LEVELS,
    COSMIC_MUTATION_CONTEXTS,
)
from hvantk.core.utils.table_utils import get_row_fields, build_rename_map, str_to_bool
from hvantk.core.utils.file_utils import resolve_compression

if TYPE_CHECKING:
    from hvantk.core.streamers.gene_catalog import GeneCatalogStreamer

logger = logging.getLogger(__name__)


def build_cosmic_cgc_submissions(
    parsed_input,
    ctx,
    **params,
):
    """Phase B builder — returns an AnnotationTable.

    Imports the COSMIC CGC TSV, renames fields, normalises tier classifications,
    optionally filters by mutation context or tier, and wraps with Provenance.

    Parameters
    ----------
    parsed_input : str | Path
        Path to the COSMIC CGC TSV file (e.g. Cosmic_Genes_v98_GRCh38.tsv.gz).
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    **params
        Optional: mutation_context (str, default "both"), min_classification (str),
                  gene_catalog (GeneCatalogStreamer), fields (list of str).
    """
    from hvantk.core.models import AnnotationTable

    mutation_context = params.get("mutation_context", "both")
    min_classification = params.get("min_classification", None)
    gene_catalog: Optional[GeneCatalogStreamer] = params.get("gene_catalog", None)
    fields = params.get("fields", None)

    if mutation_context not in COSMIC_MUTATION_CONTEXTS:
        raise ValueError(
            f"mutation_context must be one of {COSMIC_MUTATION_CONTEXTS}, "
            f"got: {mutation_context}"
        )

    if (
        min_classification is not None
        and min_classification not in COSMIC_CGC_CLASSIFICATION_LEVELS
    ):
        raise ValueError(
            f"min_classification must be one of "
            f"{COSMIC_CGC_CLASSIFICATION_LEVELS}, got: {min_classification}"
        )

    resolved_path, force_bgz = resolve_compression(str(parsed_input))

    # 1. Import
    import_kwargs = dict(
        paths=resolved_path,
        delimiter="\t",
        impute=False,
        min_partitions=10,
    )
    if force_bgz:
        import_kwargs["force_bgz"] = True
    ht = hl.import_table(**import_kwargs)

    # 2. Transform
    logger.info("Renaming COSMIC CGC fields to standardized names")
    rename_map = build_rename_map(COSMIC_CGC_FIELDS, get_row_fields(ht))
    ht = ht.rename(rename_map)

    # Normalize Tier: raw "1"/"2" -> "Tier 1"/"Tier 2"
    logger.info("Normalizing tier classification values")
    ht = ht.annotate(
        classification=hl.if_else(
            ht.classification.matches(r"^\d+$"),
            hl.literal("Tier ") + ht.classification,
            ht.classification,
        )
    )

    # Add classification_level numeric field
    classification_order = {
        level: i for i, level in enumerate(COSMIC_CGC_CLASSIFICATION_LEVELS)
    }
    ht = ht.annotate(
        classification_level=hl.literal(classification_order).get(
            ht.classification,
            hl.len(COSMIC_CGC_CLASSIFICATION_LEVELS),
        )
    )

    # Normalize boolean fields
    for bool_field in ("somatic", "germline", "hallmark"):
        if bool_field in get_row_fields(ht):
            ht = ht.annotate(**{bool_field: str_to_bool(ht[bool_field])})

    # Parse comma-separated multi-value fields into arrays
    multi_value_fields = [
        "tumour_types_somatic",
        "tumour_types_germline",
        "role_in_cancer",
        "mutation_types",
    ]
    for mv_field in multi_value_fields:
        if mv_field in get_row_fields(ht):
            ht = ht.annotate(
                **{
                    mv_field: hl.if_else(
                        hl.is_defined(ht[mv_field]) & (ht[mv_field] != ""),
                        ht[mv_field]
                        .split(",")
                        .map(lambda x: x.strip())
                        .filter(lambda x: x != ""),
                        hl.empty_array(hl.tstr),
                    )
                }
            )

    # Apply mutation_context filter
    if mutation_context == "somatic":
        logger.info("Filtering to somatic genes")
        ht = ht.filter(ht.somatic)
    elif mutation_context == "germline":
        logger.info("Filtering to germline genes")
        ht = ht.filter(ht.germline)

    # Apply min_classification filter
    if min_classification is not None:
        min_level = classification_order[min_classification]
        logger.info(
            "Filtering to classifications >= %s (level %d)",
            min_classification,
            min_level,
        )
        ht = ht.filter(ht.classification_level <= min_level)

    # Resolve gene_symbol -> hgnc_id if a gene catalog is available
    if gene_catalog is not None:
        logger.info("Resolving gene symbols to HGNC IDs via gene catalog")
        symbols = set(ht.aggregate(hl.agg.collect_as_set(ht.gene_symbol)))
        mapping = gene_catalog.map_ids(
            list(symbols), source_type="gene_symbol", target_type="hgnc_id"
        )
        mapping_literal = hl.literal(mapping)
        ht = ht.annotate(hgnc_id=mapping_literal.get(ht.gene_symbol))
        ht = ht.annotate(
            hgnc_id=hl.if_else(
                hl.is_defined(ht.hgnc_id) & ht.hgnc_id.startswith("HGNC:"),
                ht.hgnc_id.replace("HGNC:", ""),
                ht.hgnc_id,
            )
        )
        ht = ht.filter(hl.is_defined(ht.hgnc_id) & (ht.hgnc_id != ""))
        ht = ht.key_by("hgnc_id")
    else:
        logger.warning(
            "No gene catalog provided; keying by gene_symbol. "
            "Provide gene_catalog for HGNC ID resolution."
        )
        ht = ht.key_by("gene_symbol")

    # 3. Optional field selection
    if fields is not None:
        logger.info("Selecting fields: %s", fields)
        ht = ht.select(*fields)

    # 4. Wrap with provenance
    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="cosmic-cgc-v1")
    )
