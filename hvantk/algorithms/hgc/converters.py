import logging
import hail as hl
from hvantk.algorithms.hgc.constants import ADJ_GT_FIELD, VCF_EXTENSION
from hvantk.core.models.backends import algorithm, Backend

# Make gnomad import optional - only required when adjust_genotypes=True
try:
    from gnomad.utils.annotations import annotate_adj

    GNOMAD_AVAILABLE = True
except ImportError:
    GNOMAD_AVAILABLE = False
    annotate_adj = None  # defined to satisfy linters; guarded by GNOMAD_AVAILABLE


def _split_vds(
    vds: hl.vds.VariantDataset, skip_split: bool = False
) -> hl.vds.VariantDataset:
    """
    Split multi-allelic variants at the VDS level.

    This is the recommended approach for VDS data as it handles sparse
    variant/reference data correctly and produces biallelic variants.

    Parameters:
        vds: Input VariantDataset
        skip_split: If True, skip splitting

    Returns:
        Split (or unchanged) VariantDataset
    """
    if skip_split:
        logging.info("Skipping VDS-level multi-allelic split as requested.")
        return vds

    logging.info("Splitting multi-allelic variants at VDS level (sparse split)…")
    return hl.vds.split_multi(vds)


def _gt_out_of_bounds(mt: hl.MatrixTable, field: str = "GT"):
    """Predicate: this entry's call references an allele index that does not exist.

    Shared by the audit (which counts them) and the repair (which sets them missing), so the two
    can never drift apart.
    """
    gt = mt[field]
    return hl.is_defined(gt) & hl.any(
        lambda i: gt[i] >= hl.len(mt.alleles), hl.range(0, gt.ploidy)
    )


def _audit_split_variant_data(vds: hl.vds.VariantDataset) -> None:
    """Count (and report) misaligned GT/AD entries, on the SPARSE variant data.

    This is the same check that used to run on the densified MatrixTable, moved upstream. It is
    equivalent because neither defect can originate in a reference block:

    * reference-derived entries get a hom-ref call -- `to_dense_mt` either injects `hl.call(0, 0)`
      outright, or reuses a reference call that the VDS combiner already forced to be hom-ref
      (`.or_error('found reference block with non reference-genotype at' ...)`). A hom-ref call can
      never index an allele that does not exist.
    * reference-derived entries have `AD` MISSING. The combiner renames the reference block's depth
      field to `LAD`, and `to_dense_mt` fills every variant-only field with `hl.missing` on the
      reference side -- so the `hl.is_defined(AD)` guard below can never fire there.

    Both defects can therefore only live in `variant_data`, which holds just the variant records
    rather than the full samples x sites matrix. Auditing it costs a scan of the sparse data
    instead of a second full densify of the cohort.

    (Verified against Hail 0.2.137: hail/vds/methods.py `to_dense_mt`,
    hail/vds/combiner/combine.py `make_ref_entry_struct`.)

    Assumes a combiner-built VDS, which is the only kind hvantk produces. A hand-assembled VDS with
    a non-hom-ref reference block could in principle hide a defect from this count -- it would
    still be *repaired* downstream, because the repair runs on the dense matrix, but it would not
    be *reported* here.
    """
    vd = vds.variant_data
    has_gt, has_ad = "GT" in vd.entry, "AD" in vd.entry

    if not (has_gt or has_ad):
        logging.warning(
            "Skipping biallelic audit: variant data has neither 'GT' nor 'AD' "
            f"(entries: {list(vd.entry.keys())}). This is expected when the VDS-level "
            "multi-allelic split was skipped, which leaves local alleles (LGT/LAD) in place."
        )
        return

    logging.info(
        "Auditing biallelic entries on variant data (GT indices and AD lengths)…"
    )

    checks = {}
    if has_gt:
        invalid_gt = _gt_out_of_bounds(vd)
        checks["n_invalid_gt"] = hl.agg.count_where(invalid_gt)
        checks["gt_examples"] = hl.agg.filter(
            invalid_gt,
            hl.agg.take(hl.struct(locus=vd.locus, alleles=vd.alleles, GT=vd.GT), 3),
        )
    if has_ad:
        invalid_ad = hl.is_defined(vd.AD) & (hl.len(vd.AD) != hl.len(vd.alleles))
        checks["n_invalid_ad"] = hl.agg.count_where(invalid_ad)
        checks["ad_examples"] = hl.agg.filter(
            invalid_ad,
            hl.agg.take(hl.struct(locus=vd.locus, alleles=vd.alleles, AD=vd.AD), 3),
        )

    audit = vd.aggregate_entries(hl.struct(**checks))

    n_gt = audit.n_invalid_gt if has_gt else 0
    n_ad = audit.n_invalid_ad if has_ad else 0

    if n_gt == 0 and n_ad == 0:
        logging.info("✓ All entries valid (GT indices and AD lengths correct)")
        return

    if n_gt > 0:
        logging.warning(
            f"Found {n_gt} entries with out-of-bounds GT indices; they will be set to missing. "
            f"Examples: {audit.gt_examples}"
        )
    if n_ad > 0:
        # Reported only -- an AD of the wrong length is a caller/pipeline bug upstream, and
        # silently rewriting depths would hide it.
        logging.warning(
            f"Found {n_ad} entries with AD length mismatch. Examples: {audit.ad_examples}"
        )


def _apply_biallelic_gt_fix(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Set out-of-bounds genotypes to missing, as a LAZY expression on the dense matrix.

    Applied unconditionally, and that is the crux of the fix. The previous version gated this on
    `if n_invalid_gt > 0`, and *reading* that count is an eager action: it forced Hail to execute
    the whole densify plan a first time, purely to decide whether to add an annotation. Expressed
    unconditionally, the repair instead folds into the single pass that writes the MatrixTable.

    On clean data -- which is what a VDS-level split is supposed to guarantee -- this is the
    identity: `hl.if_else(False, missing, GT)` is `GT`, so the written entries are unchanged.

    It stays on the DENSE matrix on purpose: correctness of the output must not depend on the
    reference-block invariants that make the sparse audit sound.
    """
    if "GT" not in mt.entry:
        return mt

    return mt.annotate_entries(
        GT=hl.if_else(_gt_out_of_bounds(mt), hl.missing(hl.tcall), mt.GT)
    )


@algorithm(name="convert_vds_to_mt", backends=[Backend.HAIL])
def convert_vds_to_mt(
    vds_path: str,
    output_path: str,
    adjust_genotypes: bool = True,
    skip_split_multi: bool = False,
    skip_validation: bool = False,
    skip_keying_by_cols: bool = False,
    overwrite: bool = False,
) -> None:
    """
    Convert a Variant Dataset (VDS) to MatrixTable (MT).

    This function:
    1. Splits multi-allelic variants at VDS level (recommended for VDS data)
    2. Audits biallelic entries on the sparse variant data (optional, for safety)
    3. Densifies to MatrixTable
    4. Repairs out-of-bounds genotypes (lazily, fused into the write)
    5. Annotates adjusted genotypes (optional, requires gnomad)
    6. Keys by sample and writes to disk

    Parameters:
        vds_path: Path to the input VDS
        output_path: Path where the output MatrixTable will be written
        adjust_genotypes: If True, annotate with adjusted genotypes (requires gnomad)
        skip_split_multi: If True, skip splitting multi-allelic variants
        skip_validation: If True, skip both the biallelic audit and the repair
        skip_keying_by_cols: If True, skip keying the MatrixTable by columns
        overwrite: Whether to overwrite the output if it already exists

    Raises:
        RuntimeError: If adjust_genotypes=True but gnomad is not installed

    Notes:
        - VDS-level splitting is critical for correct GT/AD/PL alignment
        - After VDS split, GT/AD are already biallelic (no manual downcoding needed)
        - The densify is executed EXACTLY ONCE, by the final write. Hail is lazy and does not
          cache, so any eager action on the dense MatrixTable (an aggregate, a count) would
          re-execute the densify from the top. The audit therefore runs on the sparse variant data
          and the repair is expressed lazily; see `_audit_split_variant_data`.
        - This algorithm operates on raw `hl.MatrixTable` / `hl.VariantDataset` instances
          (genotype data). ExpressionMatrix's hail-mt backend isn't available yet (Phase J).
    """
    try:
        # Check dependencies
        if adjust_genotypes and not GNOMAD_AVAILABLE:
            raise RuntimeError(
                "adjust_genotypes=True requires the 'gnomad' package to be installed. "
                "Please install it with 'pip install gnomad' or set adjust_genotypes=False."
            )

        # Step 1: Load and split VDS
        logging.info(f"Reading VDS from {vds_path}...")
        vds = hl.vds.read_vds(vds_path)
        vds = _split_vds(vds, skip_split=skip_split_multi)

        # Step 2: Audit the SPARSE variant data, before densifying. Every eager action on the
        # dense MatrixTable re-runs the densify from scratch (Hail is lazy and does not cache), so
        # a validation aggregate there would double the cost of this whole stage.
        if skip_validation:
            logging.info("Skipping biallelic audit (skip_validation=True)")
        else:
            _audit_split_variant_data(vds)

        # Step 3: Densify to MatrixTable
        logging.info("Converting VDS to dense MatrixTable…")
        mt = hl.vds.to_dense_mt(vds)

        # Step 4: Repair out-of-bounds genotypes. A lazy expression: it costs no extra pass,
        # it fuses into the write below, and on clean data it is the identity.
        if not skip_validation:
            mt = _apply_biallelic_gt_fix(mt)

        # Step 5: Annotate adjusted genotypes (optional, requires gnomad)
        if adjust_genotypes:
            logging.info("Annotating MatrixTable with adjusted genotypes...")
            required_fields = {"GQ", "DP", "AD", "GT"}
            missing_fields = required_fields - set(mt.entry.keys())

            if missing_fields:
                logging.warning(
                    f"Cannot annotate adjusted genotypes: missing {missing_fields}. "
                    f"Available: {list(mt.entry.keys())}"
                )
            else:
                mt = annotate_adj(mt)
                logging.info("Adjusted genotype annotation completed.")

        # Step 6: Key by sample (optional)
        if not skip_keying_by_cols:
            logging.info("Keying MatrixTable by sample column 's'...")
            mt = mt.key_cols_by(mt["s"])

        # Step 7: Write output -- the ONLY action that executes the densify
        logging.info(f"Writing MatrixTable to {output_path}...")
        mt.write(output_path, overwrite=overwrite)
        logging.info("✓ MatrixTable successfully written.")

    except Exception as e:
        logging.exception("An error occurred during VDS to MT conversion.")
        raise


def convert_mt_to_multi_sample_vcf(
    mt_path: str,
    vcf_path: str,
    filter_adj_genotypes: bool = True,
    min_ac: int = 1,
    split_multi: bool = True,
) -> None:
    """
    Convert a Hail MatrixTable to a multi-sample VCF file.

    This function reads a MatrixTable, optionally filters entries to only include
    adjusted genotypes, computes variant quality control metrics, annotates the rows
    with VCF-compatible info fields, filters variants based on the minimum alternate allele
    count (AC), drops non-VCF compatible fields, and exports the result to a VCF file.

    Parameters:
        mt_path (str): Path to the input MatrixTable.
        vcf_path (str): Path where the output VCF will be written.
        filter_adj_genotypes (bool): If True, filter entries to adjusted genotypes. Recommended.
        min_ac (int): Minimum alternate allele count (AC) for a variant to be retained.
        split_multi (bool): Whether to split multi-allelic variants.
    """
    try:
        # Validate VCF path
        if vcf_path.endswith(".vcf.gz"):
            logging.warning(
                f"VCF path ends in .vcf.gz - consider using {VCF_EXTENSION} for block gzip compression. "
                f"Block gzip is the VCF standard and ensures compatibility with bcftools, tabix, GATK."
            )

        logging.info(f"Reading MatrixTable from {mt_path}...")
        mt = hl.read_matrix_table(mt_path)

        if filter_adj_genotypes:
            # check if adj field is present
            if ADJ_GT_FIELD not in mt.entry:
                logging.warning(
                    f"Cannot filter by adjusted genotypes: '{ADJ_GT_FIELD}' field not found in MatrixTable. "
                    f"This usually happens when AD field was missing during VDS→MT conversion. "
                    f"Proceeding without adj filtering. Available entry fields: {list(mt.entry.keys())}"
                )
            else:
                logging.info("Filtering entries to adjusted genotypes...")
                mt = mt.filter_entries(mt.adj, keep=True)
        else:
            logging.info("Skipping filtering for adjusted genotypes.")

        # Check current state of variants
        has_was_split = "was_split" in mt.row.keys()
        has_lpgt = "LPGT" in mt.entry.keys()

        logging.info(
            f"Pre-split check: was_split={has_was_split}, LPGT={has_lpgt}, split_multi param={split_multi}"
        )

        if not split_multi:
            logging.info("Skipping splitting multi-allelic variants: option disabled.")
        elif has_was_split:
            logging.info(
                "Skipping splitting multi-allelic variants: already split in previous step."
            )
        else:
            logging.info("Splitting multi-allelic variants now...")
            mt = hl.split_multi_hts(mt)
            logging.info("Split completed. Updating field availability check...")
            has_was_split = True
            has_lpgt = "LPGT" in mt.entry.keys()

        # CRITICAL FIX: Replace GT with LPGT BEFORE variant_qc to ensure correct computation
        # This must happen before variant_qc() because variant_qc uses GT to compute AC/AF
        if has_was_split and has_lpgt:
            logging.info(
                "Detected split variants - replacing GT with LPGT BEFORE QC for correct allele indices"
            )
            mt = mt.annotate_entries(GT=mt.LPGT)

            # Also replace AD and PL for consistency (update to local versions for split variants)
            if "LAD" in mt.entry:
                logging.info(
                    "Updating AD with LAD (local allele depths) for split variants"
                )
                mt = mt.annotate_entries(AD=mt.LAD)
            if "LPL" in mt.entry:
                logging.info(
                    "Updating PL with LPL (local phred likelihoods) for split variants"
                )
                mt = mt.annotate_entries(PL=mt.LPL)

        logging.info("Computing variant QC metrics...")
        # compute variant QC metrics and annotate/update (e.g., AC/AF) into info field
        # this recommended after filtering (e.g., adj genotypes)
        # NOW uses the correct GT (LPGT for split variants)
        mt = hl.variant_qc(mt)

        logging.info(
            "Annotating rows with VCF-compatible info fields (AF, AC, AN, call_rate)..."
        )
        mt = mt.annotate_rows(
            info=hl.struct(
                AF=mt.variant_qc.AF[1],
                AC=mt.variant_qc.AC[1],
                AN=mt.variant_qc.AN,
                call_rate=mt.variant_qc.call_rate,
            )
        )

        if min_ac >= 1:
            logging.info(f"Filtering rows with AC >= {min_ac}...")
            mt = mt.filter_rows(mt.info.AC >= min_ac, keep=True)
        else:
            logging.info("Skipping filtering rows based on AC: option disabled.")

        logging.info("Preparing MatrixTable for VCF export...")
        # Get existing entry and row fields directly from the MatrixTable
        existing_entry_fields = set(mt.entry.keys())
        existing_row_fields = set(mt.row.keys())

        # CHECKPOINT 1: Validate variant structure
        logging.info("CHECKPOINT 1: Validating variant structure before VCF export...")
        n_variants_before = mt.count_rows()

        # Check allele structure
        allele_stats = mt.aggregate_rows(
            hl.struct(
                max_alleles=hl.agg.max(hl.len(mt.alleles)),
                min_alleles=hl.agg.min(hl.len(mt.alleles)),
                n_multiallelic=hl.agg.count_where(hl.len(mt.alleles) > 2),
            )
        )
        logging.info(
            f"  Allele structure: min={allele_stats.min_alleles}, max={allele_stats.max_alleles}, "
            f"multi-allelic={allele_stats.n_multiallelic}/{n_variants_before}"
        )

        # CRITICAL: Filter to only biallelic variants
        # Even after split_multi_hts, variant_qc might have created multi-allelic entries
        if allele_stats.max_alleles > 2:
            logging.warning(
                f"Found {allele_stats.n_multiallelic} multi-allelic variants after split! Filtering to biallelic only..."
            )
            mt = mt.filter_rows(hl.len(mt.alleles) == 2)
            n_variants_after = mt.count_rows()
            logging.info(
                f"  Filtered: {n_variants_before} → {n_variants_after} variants (removed {n_variants_before - n_variants_after})"
            )

        # CHECKPOINT 2: Verify split variant field replacements
        logging.info("CHECKPOINT 2: Verifying split variant field replacements...")
        if "was_split" in existing_row_fields:
            logging.info("  Split variants detected (was_split flag present)")
            # Check if GT has been replaced with LPGT (should have happened before variant_qc)
            if "LPGT" in existing_entry_fields:
                logging.info(
                    "  ✓ GT was replaced with LPGT before variant_qc (correct)"
                )
            else:
                logging.warning(
                    "  WARNING: was_split present but LPGT not found - this may cause issues"
                )
        else:
            logging.info("  No split variants detected (was_split flag not present)")

        # CHECKPOINT 3: Validate genotype indices
        logging.info("CHECKPOINT 3: Validating genotype allele indices...")
        # Check if any GT has alleles > 1 (should only be 0 or 1 for biallelic)
        # Use hl.call.unphased_diploid_gt_index_call to properly extract allele indices
        gt_validation = mt.aggregate_entries(
            hl.struct(
                n_defined=hl.agg.count_where(hl.is_defined(mt.GT)),
                n_invalid=hl.agg.count_where(
                    hl.is_defined(mt.GT)
                    & (
                        (
                            mt.GT.unphased_diploid_gt_index() >= 3
                        )  # For biallelic: 0/0=0, 0/1=1, 1/1=2, anything >=3 is invalid
                    )
                ),
                example_invalid=hl.agg.filter(
                    hl.is_defined(mt.GT) & (mt.GT.unphased_diploid_gt_index() >= 3),
                    hl.agg.take(
                        hl.struct(locus=mt.locus, alleles=mt.alleles, GT=mt.GT), 5
                    ),
                ),
            )
        )

        logging.info(
            f"  GT validation: {gt_validation.n_defined} defined genotypes, {gt_validation.n_invalid} invalid"
        )

        if gt_validation.n_invalid > 0:
            logging.error(
                f"  ERROR: Found {gt_validation.n_invalid} genotypes with invalid allele indices!"
            )
            logging.error(
                f"  Example invalid genotypes: {gt_validation.example_invalid}"
            )
            logging.error(
                "  These genotypes reference allele indices that don't exist in biallelic variants"
            )
            logging.error("  Setting these genotypes to missing...")

            # Set invalid genotypes to missing
            # A valid biallelic genotype should have GT index 0 (0/0), 1 (0/1), or 2 (1/1)
            mt = mt.annotate_entries(
                GT=hl.if_else(
                    hl.is_defined(mt.GT) & (mt.GT.unphased_diploid_gt_index() >= 3),
                    hl.missing(hl.tcall),
                    mt.GT,
                )
            )

            # Verify the fix
            remaining_invalid = mt.aggregate_entries(
                hl.agg.count_where(
                    hl.is_defined(mt.GT) & (mt.GT.unphased_diploid_gt_index() >= 3)
                )
            )
            logging.info(
                f"  After filtering: {remaining_invalid} invalid genotypes remaining"
            )
        else:
            logging.info("  ✓ All genotypes have valid allele indices")

        # CHECKPOINT 4: Select VCF-compatible fields
        logging.info("CHECKPOINT 4: Selecting VCF-compatible entry fields...")
        vcf_standard_entry_fields = {"GT", "DP", "GQ", "PID", "SB"}

        # Add AD and PL if they exist in the MT
        if "AD" in mt.entry.keys():
            vcf_standard_entry_fields.add("AD")
        if "PL" in mt.entry.keys():
            vcf_standard_entry_fields.add("PL")

        entry_fields_to_keep = [
            f for f in vcf_standard_entry_fields if f in mt.entry.keys()
        ]

        if entry_fields_to_keep:
            logging.info(
                f"  Keeping entry fields: {', '.join(sorted(entry_fields_to_keep))}"
            )
            mt = mt.select_entries(*entry_fields_to_keep)
        else:
            logging.warning("  No standard VCF entry fields found")

        # CHECKPOINT 5: Drop incompatible row fields
        logging.info("CHECKPOINT 5: Dropping VCF-incompatible row fields...")
        row_fields_to_drop = {"variant_qc", "a_index", "was_split"}
        row_fields_present = [f for f in row_fields_to_drop if f in mt.row.keys()]
        if row_fields_present:
            logging.info(f"  Dropping row fields: {', '.join(row_fields_present)}")
            mt = mt.drop(*row_fields_present)

        # CHECKPOINT 6: Final validation before export
        logging.info("CHECKPOINT 6: Final validation before VCF export...")
        final_counts = mt.count()
        logging.info(
            f"  Final MatrixTable: {final_counts[0]} variants × {final_counts[1]} samples"
        )

        # Validate final state
        final_stats = mt.aggregate_rows(
            hl.struct(
                n_biallelic=hl.agg.count_where(hl.len(mt.alleles) == 2),
                n_multiallelic=hl.agg.count_where(hl.len(mt.alleles) > 2),
                total=hl.agg.count(),
            )
        )
        logging.info(
            f"  Variants: biallelic={final_stats.n_biallelic}, multi-allelic={final_stats.n_multiallelic}"
        )

        if final_stats.n_multiallelic > 0:
            raise ValueError(
                f"Cannot export VCF: still have {final_stats.n_multiallelic} multi-allelic variants! "
                "VCF export requires all variants to be biallelic."
            )

        logging.info("  ✓ All validations passed")

        logging.info(f"Exporting VCF to {vcf_path}...")
        hl.export_vcf(mt, vcf_path)
        logging.info("VCF successfully written.")

    except Exception as e:
        logging.exception("An error occurred during conversion to VCF.")
        raise
