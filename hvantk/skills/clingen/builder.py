"""Hail Table builder for the ClinGen Gene-Disease Validity resource.

Owns the Phase B ``build_clingen_gene_disease`` builder. Turns the ClinGen
Gene-Disease Validity CSV into an ``AnnotationTable`` keyed by
``(hgnc_id, mondo_id)``. ClinGen provides curated gene-disease associations
with evidence-based classifications (Definitive, Strong, Moderate, Limited,
Disputed, Refuted, No Known Disease Relationship).

ClinGen CSV files have a varying metadata header before the actual data
header (the line starting with ``GENE SYMBOL``) and ``++++++`` separator
rows; the builder preprocesses these out before passing the file to
``hl.import_table``.

``cleanup_temp_file`` lives in ``hvantk.core.utils.hail_helpers`` as shared
infrastructure.
"""

from __future__ import annotations

import logging

import hail as hl

from hvantk.skills.clingen.shared.constants import (
    CLINGEN_CLASSIFICATION_LEVELS,
    CLINGEN_GENE_DISEASE_FIELDS,
)
from hvantk.core.utils.hail_helpers import cleanup_temp_file
from hvantk.core.utils.table_utils import (
    get_row_fields,
    strip_curie_prefix,
    annotate_classification_level,
    filter_min_classification,
)

logger = logging.getLogger(__name__)


def build_clingen_gene_disease(
    parsed_input,
    ctx,
    *,
    min_classification=None,
    fields=None,
):
    """Phase B builder — returns an AnnotationTable keyed by
    ``(hgnc_id, mondo_id)``.
    """
    from hvantk.core.models import AnnotationTable

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
            ht = ht.annotate(hgnc_id=strip_curie_prefix(ht.hgnc_id, "HGNC:"))
        if "mondo_id" in row_fields:
            ht = ht.annotate(mondo_id=strip_curie_prefix(ht.mondo_id, "MONDO:"))

        ht = annotate_classification_level(ht, CLINGEN_CLASSIFICATION_LEVELS)
        if min_classification is not None:
            ht = filter_min_classification(
                ht, CLINGEN_CLASSIFICATION_LEVELS, min_classification
            )

        ht = ht.key_by("hgnc_id", "mondo_id")

        if fields is not None:
            ht = ht.select(*fields)

        return AnnotationTable.from_hail(
            ht, provenance=ctx.provenance(schema_id="clingen-gene-disease-v1")
        )
    except Exception:
        cleanup_temp_file(tmp_path)
        raise
    # NOTE: tmp_path is intentionally NOT cleaned up on success — Hail's lazy
    # evaluation may read from it later when artifact.save() materializes the table.
    # The OS will clean up the Hail temp directory at session end.
