"""Hail Table builder for the dbNSFP variant functional annotation database.

Owns the Phase B ``build_dbnsfp_variants`` builder. Imports the dbNSFP
TSV/BGZ, parses variant coordinates to ``(locus, alleles)``, optionally
groups transcript scores and common prefixes into structs, and wraps with
Provenance.
"""
from __future__ import annotations

import logging

import hail as hl

from hvantk.core.utils.table_utils import get_row_fields
from hvantk.core.utils.file_utils import resolve_compression

logger = logging.getLogger(__name__)


def build_dbnsfp_variants(
    parsed_input,
    ctx,
    **params,
):
    """Phase B builder — returns an AnnotationTable.

    Imports the dbNSFP TSV/BGZ, parses variant coordinates to (locus, alleles),
    optionally groups transcript scores and common prefixes into structs, and
    wraps with Provenance.

    Parameters
    ----------
    parsed_input : str | Path
        Path to the dbNSFP TSV/BGZ input file.
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    **params
        Optional: reference_genome (str, default "GRCh38"),
                  min_partitions (int, default 200),
                  force_bgz (bool, default True),
                  parse_transcript_scores (bool, default True),
                  group_prefixes (list of str),
                  auto_convert_bgz (bool, default False).
    """
    from hvantk.core.models import AnnotationTable

    reference_genome = params.get("reference_genome", "GRCh38")
    min_partitions = params.get("min_partitions", 200)
    force_bgz = params.get("force_bgz", True)
    parse_transcript_scores = params.get("parse_transcript_scores", True)
    group_prefixes = params.get("group_prefixes", None)
    auto_convert_bgz = params.get("auto_convert_bgz", False)

    # Resolve compression: detect gz vs bgzf, optionally convert
    input_path_resolved, force_bgz = resolve_compression(
        str(parsed_input),
        force_bgz=force_bgz,
        auto_convert=auto_convert_bgz,
    )

    ht = hl.import_table(
        paths=input_path_resolved,
        min_partitions=min_partitions,
        impute=False,
        missing=".",
        force_bgz=force_bgz,
    )

    # Normalize chromosome field and construct variant key
    row_fields = get_row_fields(ht)
    if "#chr" in row_fields:
        ht = ht.rename({"#chr": "chr"})
    else:
        if "chr" not in row_fields:
            raise ValueError("dbNSFP input missing '#chr' or 'chr' column")

    _chr_str = hl.str(ht["chr"])
    ht = ht.annotate(
        chr=hl.if_else(
            _chr_str.lower().startswith("chr"), _chr_str, hl.str("chr") + _chr_str
        )
    )

    # Build variant_key: chr:pos:ref:alt
    row_fields = get_row_fields(ht)
    if (
        "pos(1-based)" not in row_fields
        or "ref" not in row_fields
        or "alt" not in row_fields
    ):
        raise ValueError(
            "dbNSFP input missing required columns: 'pos(1-based)', 'ref', or 'alt'"
        )

    variant_key_expr = hl.array(
        [ht.chr, hl.str(ht["pos(1-based)"]), ht.ref, ht.alt]
    )
    ht = ht.annotate(variant_key=hl.delimit(variant_key_expr, ":"))

    # Parse to locus/alleles
    ht = ht.annotate(
        **hl.parse_variant(ht.variant_key, reference_genome=reference_genome)
    )

    # Key the table and cleanup staging columns
    ht = ht.key_by("locus", "alleles")
    ht = ht.drop("variant_key", "chr", "pos(1-based)", "ref", "alt")

    # Transcript-specific score parsing
    row_fields = get_row_fields(ht)
    if parse_transcript_scores and "Ensembl_transcriptid" in row_fields:
        logger.info(
            "Parsing transcript-specific scores into dicts keyed by Ensembl_transcriptid"
        )
        ht = ht.annotate(Ensembl_transcriptid=hl.str(ht.Ensembl_transcriptid))
        ht = ht.annotate(Ensembl_transcriptid=ht.Ensembl_transcriptid.split(";"))

        row_fields_list = list(get_row_fields(ht))
        score_fields = [
            f for f in row_fields_list if f.endswith("_score") or f == "CADD_phred"
        ]

        def _to_float_array(s):
            s_def = hl.or_else(s, "")
            arr = s_def.split(";")
            return hl.map(hl.parse_float, arr)

        def _single_to_dict(val):
            return hl.dict(
                hl.zip(
                    ht.Ensembl_transcriptid,
                    hl.map(lambda _x: hl.parse_float(val), ht.Ensembl_transcriptid),
                )
            )

        ann = {}
        for f in score_fields:
            is_multi = hl.is_defined(ht[f]) & ht[f].contains(";")
            ann[f] = hl.if_else(
                is_multi,
                hl.dict(hl.zip(ht.Ensembl_transcriptid, _to_float_array(ht[f]))),
                _single_to_dict(ht[f]),
            )
        if ann:
            ht = ht.annotate(**ann)

    # Group common prefixes into structs
    prefixes = group_prefixes or [
        "gnomAD",
        "ExAC",
        "1000Gp3",
        "ESP6500",
        "clinvar",
    ]
    for prefix in prefixes:
        row_fields_list = list(get_row_fields(ht))
        pref_fields = [
            f for f in row_fields_list if f != prefix and f.startswith(prefix)
        ]
        if pref_fields:
            logger.info(
                "Grouping %s* fields into struct '%s' (%d fields)",
                prefix,
                prefix,
                len(pref_fields),
            )
            ht = ht.annotate(
                **{prefix: hl.struct(**{f: ht[f] for f in pref_fields})}
            )
            ht = ht.drop(*pref_fields)

    # 3. Wrap with provenance
    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="dbnsfp-v1")
    )
