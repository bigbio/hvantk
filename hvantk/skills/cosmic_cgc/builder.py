"""Hail Table builder for the COSMIC Cancer Gene Census (CGC) resource.

Plugin builder: imports the COSMIC CGC TSV, normalises tier classifications,
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
    COSMIC_CGC_HEADER_GENERATIONS,
    COSMIC_MUTATION_CONTEXTS,
)
from hvantk.core.utils.table_utils import (
    get_row_fields,
    build_rename_map,
    str_to_bool,
    annotate_classification_level,
    filter_min_classification,
)
from hvantk.core.utils.file_utils import resolve_compression

if TYPE_CHECKING:
    from hvantk.core.streamers.gene_catalog import GeneCatalogStreamer

logger = logging.getLogger(__name__)


def _header_generation(row_fields: set[str]) -> str:
    """Return the schema ID of the header generation a renamed table carries.

    Picks the generation in ``COSMIC_CGC_HEADER_GENERATIONS`` with the fewest
    missing columns, then the fewest columns it does not define. Raises
    ValueError if that generation still lacks a column, so a partly renamed
    header is never stamped with a schema ID whose fields it does not have.
    """

    def distance(schema_id: str) -> tuple[int, int]:
        fields = set(COSMIC_CGC_HEADER_GENERATIONS[schema_id].values())
        return len(fields - row_fields), len(row_fields - fields)

    schema_id = min(COSMIC_CGC_HEADER_GENERATIONS, key=distance)
    field_map = COSMIC_CGC_HEADER_GENERATIONS[schema_id]
    missing = [raw for raw, field in field_map.items() if field not in row_fields]
    extra = sorted(row_fields - set(field_map.values()))
    if missing:
        raise ValueError(
            f"COSMIC CGC: the header is closest to the {schema_id} generation but "
            f"lacks {', '.join(missing)} (columns it does not define: "
            f"{', '.join(extra) or 'none'}). A renamed column or a new header "
            "generation needs its own column map; see "
            "hvantk/skills/cosmic_cgc/SKILL.md s8."
        )
    if extra:
        logger.warning(
            "COSMIC CGC: keeping column(s) the %s header does not define: %s",
            schema_id,
            ", ".join(extra),
        )
    return schema_id


def build_cosmic_cgc_submissions(
    parsed_input,
    ctx,
    **params,
):
    """Plugin builder — returns an AnnotationTable.

    Imports the COSMIC CGC TSV, renames fields, normalises tier classifications,
    optionally filters by mutation context or tier, and wraps with Provenance.

    Parameters
    ----------
    parsed_input : str | Path
        Path to the COSMIC CGC TSV file (e.g.
        Cosmic_CancerGeneCensus_v10x_GRCh38.tsv.gz).
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    **params
        Optional: mutation_context (str, default "both"), min_classification (str),
                  gene_catalog (GeneCatalogStreamer), fields (list of str).

    Raises
    ------
    ValueError
        If the header lacks a column of the generation it is closest to, or if
        gene_catalog resolves none of the table's gene symbols.
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
    else:
        # The real COSMIC CGC export is standard (non-block) gzip, not BGZF
        # (confirmed from its magic bytes). hl.import_table refuses to read
        # ANY non-BGZF .gz path at all without this flag -- regardless of
        # min_partitions (confirmed empirically with it unset, 1, and 10) --
        # raising HailException: ".gz cannot be loaded in parallel". Harmless
        # for uncompressed input: partition count is unaffected either way,
        # and resolve_compression() already separates genuine bgzf (handled
        # above) from plain gzip and uncompressed, both of which land here.
        import_kwargs["force"] = True
    ht = hl.import_table(**import_kwargs)

    # 2. Transform
    logger.info("Renaming COSMIC CGC fields to standardized names")
    rename_map = build_rename_map(COSMIC_CGC_FIELDS, get_row_fields(ht))
    ht = ht.rename(rename_map)
    schema_id = _header_generation(get_row_fields(ht))

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
    ht = annotate_classification_level(ht, COSMIC_CGC_CLASSIFICATION_LEVELS)

    # Normalize boolean fields
    for bool_field in ("somatic", "germline", "hallmark"):
        if bool_field in get_row_fields(ht):
            ht = ht.annotate(**{bool_field: str_to_bool(ht[bool_field])})

    # Cast genome coordinate fields. Only the current (v103+) header has them,
    # and _header_generation() requires both there, so for a current export
    # the cast and the count of non-integer values always run.
    # hl.parse_int32 is missing-tolerant: an empty or non-numeric value
    # becomes missing rather than raising (confirmed: 6/763 rows in a
    # licensed v103 export leave these fields empty).
    for coord_field in ("genome_start", "genome_stop"):
        if coord_field in get_row_fields(ht):
            raw_coordinate = ht[coord_field]
            parsed_coordinate = hl.parse_int32(raw_coordinate)
            invalid_count = ht.aggregate(
                hl.agg.count_where(
                    hl.is_defined(raw_coordinate)
                    & (raw_coordinate != "")
                    & hl.is_missing(parsed_coordinate)
                )
            )
            if invalid_count:
                logger.warning(
                    "COSMIC CGC: %d non-empty %s value(s) could not be parsed as int32",
                    invalid_count,
                    coord_field,
                )
            ht = ht.annotate(**{coord_field: parsed_coordinate})

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
        logger.info("Filtering to classifications >= %s", min_classification)
        ht = filter_min_classification(
            ht, COSMIC_CGC_CLASSIFICATION_LEVELS, min_classification
        )

    # Resolve gene_symbol -> hgnc_id if a gene catalog is available
    if gene_catalog is not None:
        logger.info("Resolving gene symbols to HGNC IDs via gene catalog")
        # collect_as_set keeps missing and empty values; neither is a symbol.
        symbols = set(ht.aggregate(hl.agg.collect_as_set(ht.gene_symbol))) - {None, ""}
        mapping = gene_catalog.map_ids(
            list(symbols), source_type="gene_symbol", target_type="hgnc_id"
        )
        mapping = {symbol: hgnc_id for symbol, hgnc_id in mapping.items() if hgnc_id}
        unresolved = sorted(symbols - mapping.keys())
        if symbols and len(unresolved) == len(symbols):
            raise ValueError(
                f"COSMIC CGC: the gene catalog resolved none of the {len(symbols)} "
                f"gene symbols to an HGNC ID (e.g. {', '.join(unresolved[:5])}); "
                "check that it is an HGNC lookup table"
            )
        if unresolved:
            logger.warning(
                "COSMIC CGC: resolved %d of %d gene symbols to HGNC IDs; dropping "
                "the rows of %d unresolved symbol(s), e.g. %s",
                len(symbols) - len(unresolved),
                len(symbols),
                len(unresolved),
                ", ".join(unresolved[:10]),
            )
        else:
            logger.info(
                "COSMIC CGC: resolved all %d gene symbols to HGNC IDs", len(symbols)
            )
        mapping_literal = hl.literal(mapping, hl.tdict(hl.tstr, hl.tstr))
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

    # 4. Wrap with provenance (schema ID of the detected header generation)
    return AnnotationTable.from_hail(ht, provenance=ctx.provenance(schema_id=schema_id))
