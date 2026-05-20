"""Hail Table builder for the ClinGen Gene-Disease Validity resource.

This module owns ``create_clingen_gene_disease_tb``, the canonical builder
that turns the ClinGen Gene-Disease Validity CSV into a Hail Table. It was
migrated out of :mod:`hvantk.core.builders.table` so that everything
ClinGen-specific (builder, downloader, dataset class, tests, fixtures,
SKILL) lives under the plugin folder at :mod:`hvantk.skills.clingen`.

The shared helpers ``_create_table_base`` and ``_cleanup_temp_file`` and the
``get_row_fields`` utility intentionally stay in their existing modules
because they are reused by other builders.
"""

from __future__ import annotations

import logging
from typing import List, Optional

import hail as hl

from hvantk.core.constants import (
    CLINGEN_CLASSIFICATION_LEVELS,
    CLINGEN_GENE_DISEASE_FIELDS,
)
from hvantk.core.builders.table import _cleanup_temp_file, _create_table_base
from hvantk.core.utils.table_utils import get_row_fields

logger = logging.getLogger(__name__)


def create_clingen_gene_disease_tb(
    input_path: str,
    output_path: str,
    key_by: str = "gene_disease",
    min_classification: Optional[str] = None,
    fields: Optional[List[str]] = None,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> "hl.Table":
    """
    Create a Hail Table from a ClinGen Gene-Disease Validity CSV file.

    ClinGen provides curated gene-disease associations with evidence-based
    classifications (Definitive, Strong, Moderate, Limited, etc.).

    Example usage:
        # Default: keyed by (hgnc_id, mondo_id)
        ht = create_clingen_gene_disease_tb(
            input_path="/path/to/clingen.csv",
            output_path="/path/to/output.ht"
        )

        # Gene-level aggregation (for joining with other gene tables)
        ht = create_clingen_gene_disease_tb(
            input_path="/path/to/clingen.csv",
            output_path="/path/to/output.ht",
            key_by="gene",
            min_classification="Moderate"
        )

    Parameters
    ----------
    input_path : str
        Path to the ClinGen Gene-Disease Validity CSV file.
    output_path : str
        Path to write the output Hail Table.
    key_by : str, optional
        Keying strategy:
        - "gene_disease" (default): Key by (hgnc_id, mondo_id) - preserves full granularity
        - "gene": Aggregate diseases per gene, key by hgnc_id
    min_classification : str, optional
        Filter to classifications at or above this level. Valid values:
        "Definitive", "Strong", "Moderate", "Limited", "Disputed", "Refuted".
        If None, includes all classifications (default: None).
    fields : list of str, optional
        List of fields to select from the table (default: None, keeps all).
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).

    Returns
    -------
    hl.Table
        Hail Table with ClinGen gene-disease validity annotations.

    Notes
    -----
    The ClinGen CSV file has a 6-line metadata header that is skipped during import.

    Classification hierarchy (from strongest to weakest):
    1. Definitive
    2. Strong
    3. Moderate
    4. Limited
    5. Disputed
    6. Refuted
    7. No Known Disease Relationship

    When key_by="gene", diseases are aggregated per gene with fields:
    - disease_labels: set of disease labels
    - disease_mondo_pairs: set of (disease_label, mondo_id) pairs
    - mondo_ids: set of MONDO IDs
    - classifications: set of classification levels
    - max_classification_level: numeric level of highest classification
    - max_classification_label: label of highest classification
    - n_diseases: count of associated diseases
    """
    if key_by not in ("gene_disease", "gene"):
        raise ValueError(f"key_by must be 'gene_disease' or 'gene', got: {key_by}")

    if (
        min_classification is not None
        and min_classification not in CLINGEN_CLASSIFICATION_LEVELS
    ):
        raise ValueError(
            f"min_classification must be one of {CLINGEN_CLASSIFICATION_LEVELS}, "
            f"got: {min_classification}"
        )

    # Preprocess CSV to extract header and data rows
    # ClinGen files have varying metadata lines before the actual header
    # The header row starts with "GENE SYMBOL" and separator rows contain "++++++"
    # Use Hadoop API to support cloud URIs (gs://, s3://) and distributed Spark clusters
    logger.info("Preprocessing ClinGen CSV: finding header and filtering metadata")

    # Use Hail's temp file utility to create a Hadoop-accessible temp file
    tmp_path = hl.utils.new_temp_file(prefix="clingen_", extension=".csv")
    try:
        # Stream the file line-by-line to find header and skip metadata/separators
        with hl.hadoop_open(input_path, "r") as f:
            with hl.hadoop_open(tmp_path, "w") as out:
                found_header = False
                for line in f:
                    # Skip separator lines (contain "++++++")
                    if "++++++" in line:
                        continue
                    # Look for header row (starts with "GENE SYMBOL" in quotes or unquoted)
                    if not found_header:
                        if '"GENE SYMBOL"' in line or line.startswith("GENE SYMBOL"):
                            found_header = True
                            out.write(line)
                        # Skip metadata lines before header
                        continue
                    # Write all data lines after header
                    out.write(line)

                if not found_header:
                    raise RuntimeError(
                        f'ClinGen header "GENE SYMBOL" not found in {input_path}'
                    )

        logger.info(f"Preprocessed file written to {tmp_path}")

    except Exception as e:
        # Clean up temp file if preprocessing fails
        _cleanup_temp_file(tmp_path)
        raise RuntimeError(f"Failed to preprocess ClinGen CSV: {e}") from e

    try:

        def import_func():
            return hl.import_table(
                paths=tmp_path,
                delimiter=",",
                quote='"',
                impute=False,
                min_partitions=10,
            )

        def transform(ht: hl.Table) -> hl.Table:
            # Rename fields to standardized names
            logger.info("Renaming fields to standardized names")
            rename_map = {
                k: v
                for k, v in CLINGEN_GENE_DISEASE_FIELDS.items()
                if k in get_row_fields(ht)
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
                level: i for i, level in enumerate(CLINGEN_CLASSIFICATION_LEVELS)
            }
            ht = ht.annotate(
                classification_level=hl.literal(classification_order).get(
                    ht.classification, hl.len(CLINGEN_CLASSIFICATION_LEVELS)
                )
            )

            # Apply min_classification filter if specified
            if min_classification is not None:
                min_level = classification_order[min_classification]
                logger.info(
                    f"Filtering to classifications >= {min_classification} (level {min_level})"
                )
                ht = ht.filter(ht.classification_level <= min_level)

            # Apply keying strategy
            if key_by == "gene_disease":
                logger.info("Keying by (hgnc_id, mondo_id)")
                ht = ht.key_by("hgnc_id", "mondo_id")
            else:  # key_by == "gene"
                logger.info("Aggregating by gene (hgnc_id)")
                ht = (
                    ht.group_by("hgnc_id", "gene_symbol")
                    .aggregate(
                        disease_labels=hl.agg.collect_as_set(ht.disease_label),
                        disease_mondo_pairs=hl.agg.collect_as_set(
                            hl.struct(
                                disease_label=ht.disease_label, mondo_id=ht.mondo_id
                            )
                        ),
                        mondo_ids=hl.agg.collect_as_set(ht.mondo_id),
                        classifications=hl.agg.collect_as_set(ht.classification),
                        modes_of_inheritance=hl.agg.collect_as_set(
                            ht.mode_of_inheritance
                        ),
                        max_classification_level=hl.agg.min(ht.classification_level),
                        n_diseases=hl.agg.count(),
                    )
                    .key_by("hgnc_id")
                )
                # Add max classification label with safe index clamping
                # max_classification_level can be len(CLINGEN_CLASSIFICATION_LEVELS) for
                # unknown classifications, so we clamp to valid index range
                classification_labels = hl.literal(CLINGEN_CLASSIFICATION_LEVELS)
                safe_index = hl.min(
                    ht.max_classification_level, hl.len(classification_labels) - 1
                )
                ht = ht.annotate(
                    max_classification_label=classification_labels[safe_index]
                )

            return ht

        clingen_tb = _create_table_base(
            source_name="ClinGen Gene-Disease Validity",
            input_path=input_path,
            output_path=output_path,
            import_func=import_func,
            transform_func=transform,
            fields=fields,
            overwrite=overwrite,
            export_tsv=export_tsv,
        )

        return clingen_tb

    finally:
        # Clean up temp file
        _cleanup_temp_file(tmp_path)


def build_clingen_gene_disease(
    parsed_input,
    ctx,
    *,
    key_by: str = "gene_disease",
    min_classification=None,
    fields=None,
):
    """Phase B builder — returns an AnnotationTable.

    Preserves the CSV header preprocessing step from create_clingen_gene_disease_tb
    (ClinGen files have a varying metadata header + ++++++ separator lines).
    Only supports key_by='gene_disease' for now; the 'gene' aggregation mode
    still requires the legacy function.
    """
    from hvantk.core.models import AnnotationTable

    if key_by != "gene_disease":
        raise NotImplementedError(
            f"Phase B builder currently supports key_by='gene_disease' only; "
            f"got {key_by}. Use the legacy create_clingen_gene_disease_tb for now."
        )
    if (
        min_classification is not None
        and min_classification not in CLINGEN_CLASSIFICATION_LEVELS
    ):
        raise ValueError(
            f"min_classification must be one of {CLINGEN_CLASSIFICATION_LEVELS}, "
            f"got: {min_classification}"
        )

    # Preprocess CSV (skip metadata + separator rows; isolate the GENE SYMBOL header)
    tmp_path = hl.utils.new_temp_file(prefix="clingen_", extension=".csv")
    try:
        with hl.hadoop_open(str(parsed_input), "r") as f:
            with hl.hadoop_open(tmp_path, "w") as out:
                found_header = False
                for line in f:
                    if "++++++" in line:
                        continue
                    if not found_header:
                        if '"GENE SYMBOL"' in line or line.startswith("GENE SYMBOL"):
                            found_header = True
                            out.write(line)
                        continue
                    out.write(line)
                if not found_header:
                    raise RuntimeError(
                        f'ClinGen header "GENE SYMBOL" not found in {parsed_input}'
                    )

        ht = hl.import_table(
            paths=tmp_path,
            delimiter=",",
            quote='"',
            impute=False,
            min_partitions=10,
        )

        rename_map = {
            k: v
            for k, v in CLINGEN_GENE_DISEASE_FIELDS.items()
            if k in get_row_fields(ht)
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
            level: i for i, level in enumerate(CLINGEN_CLASSIFICATION_LEVELS)
        }
        ht = ht.annotate(
            classification_level=hl.literal(classification_order).get(
                ht.classification, hl.len(CLINGEN_CLASSIFICATION_LEVELS)
            )
        )

        if min_classification is not None:
            min_level = classification_order[min_classification]
            ht = ht.filter(ht.classification_level <= min_level)

        ht = ht.key_by("hgnc_id", "mondo_id")

        if fields is not None:
            ht = ht.select(*fields)

        return AnnotationTable.from_hail(
            ht, provenance=ctx.provenance(schema_id="clingen-gene-disease-v1")
        )
    except Exception:
        _cleanup_temp_file(tmp_path)
        raise
    # NOTE: tmp_path is intentionally NOT cleaned up on success — Hail's lazy
    # evaluation may read from it later when artifact.save() materializes the table.
    # The OS will clean up the Hail temp directory at session end.
