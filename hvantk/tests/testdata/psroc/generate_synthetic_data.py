#!/usr/bin/env python
"""
Generate synthetic test data for PSROC pipeline testing.

This script creates reproducible synthetic ClinVar-like and dbNSFP-like
datasets suitable for running the complete PSROC pipeline including
ROC curve generation and plot visualization.

Usage:
    python generate_synthetic_data.py [--output-dir DIR] [--seed SEED]

The generated data includes:
- 100 variants across BRCA1, BRCA2, and TP53 genes
- 50 pathogenic, 40 benign, and 10 VUS variants
- 4 prediction scores with varying discriminative power and missingness
"""

import argparse
import csv
import random
from pathlib import Path
from typing import Optional


# Gene coordinates (GRCh38)
GENE_INFO = {
    "BRCA1": {"chr": "chr17", "start": 43044295, "end": 43170245},
    "BRCA2": {"chr": "chr13", "start": 32315086, "end": 32400266},
    "TP53": {"chr": "chr17", "start": 7668402, "end": 7687538},
}

# ClinVar significance labels
PATHOGENIC_LABELS = ["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"]
BENIGN_LABELS = ["Benign", "Likely_benign", "Benign/Likely_benign"]
VUS_LABELS = ["Uncertain_significance"]

# Review status values (maps to star ratings)
REVIEW_STATUS_1_STAR = ["criteria_provided,_single_submitter"]
REVIEW_STATUS_2_STAR = ["criteria_provided,_multiple_submitters,_no_conflicts"]


def generate_variant_position(gene: str, used_positions: set, rng: random.Random) -> int:
    """Generate a unique position within a gene's genomic range."""
    info = GENE_INFO[gene]
    for _ in range(1000):  # Prevent infinite loop
        pos = rng.randint(info["start"], info["end"])
        if pos not in used_positions:
            used_positions.add(pos)
            return pos
    raise ValueError(f"Could not generate unique position for {gene}")


def generate_ref_alt(rng: random.Random) -> tuple[str, str]:
    """Generate random ref/alt alleles (SNVs only for simplicity)."""
    bases = ["A", "C", "G", "T"]
    ref = rng.choice(bases)
    alt = rng.choice([b for b in bases if b != ref])
    return ref, alt


def generate_score(
    is_pathogenic: bool,
    path_mean: float,
    path_std: float,
    benign_mean: float,
    benign_std: float,
    min_val: float,
    max_val: float,
    rng: random.Random,
) -> float:
    """Generate a prediction score with class-specific distribution."""
    if is_pathogenic:
        value = rng.gauss(path_mean, path_std)
    else:
        value = rng.gauss(benign_mean, benign_std)
    return max(min_val, min(max_val, value))


def generate_synthetic_clinvar(
    output_path: str,
    seed: int = 42,
    n_pathogenic: int = 50,
    n_benign: int = 40,
    n_vus: int = 10,
) -> list[dict]:
    """
    Generate synthetic ClinVar-like TSV data.

    Args:
        output_path: Path to write the TSV file
        seed: Random seed for reproducibility
        n_pathogenic: Number of pathogenic variants
        n_benign: Number of benign variants
        n_vus: Number of VUS variants (will be excluded in PSROC)

    Returns:
        List of variant dictionaries for use in dbNSFP generation
    """
    rng = random.Random(seed)
    variants = []
    used_positions: set[int] = set()

    genes = list(GENE_INFO.keys())

    # Generate pathogenic variants
    for i in range(n_pathogenic):
        gene = genes[i % len(genes)]
        pos = generate_variant_position(gene, used_positions, rng)
        ref, alt = generate_ref_alt(rng)
        clnsig = rng.choice(PATHOGENIC_LABELS)
        review = rng.choice(REVIEW_STATUS_1_STAR + REVIEW_STATUS_2_STAR)

        variants.append({
            "chr": GENE_INFO[gene]["chr"],
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "CLNSIG": clnsig,
            "CLNREVSTAT": review,
            "GENEINFO": f"{gene}:1234",
            "is_pathogenic": True,
            "gene": gene,
        })

    # Generate benign variants
    for i in range(n_benign):
        gene = genes[i % len(genes)]
        pos = generate_variant_position(gene, used_positions, rng)
        ref, alt = generate_ref_alt(rng)
        clnsig = rng.choice(BENIGN_LABELS)
        review = rng.choice(REVIEW_STATUS_1_STAR + REVIEW_STATUS_2_STAR)

        variants.append({
            "chr": GENE_INFO[gene]["chr"],
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "CLNSIG": clnsig,
            "CLNREVSTAT": review,
            "GENEINFO": f"{gene}:1234",
            "is_pathogenic": False,
            "gene": gene,
        })

    # Generate VUS variants (will be excluded)
    for i in range(n_vus):
        gene = genes[i % len(genes)]
        pos = generate_variant_position(gene, used_positions, rng)
        ref, alt = generate_ref_alt(rng)
        clnsig = rng.choice(VUS_LABELS)
        review = rng.choice(REVIEW_STATUS_1_STAR)

        variants.append({
            "chr": GENE_INFO[gene]["chr"],
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "CLNSIG": clnsig,
            "CLNREVSTAT": review,
            "GENEINFO": f"{gene}:1234",
            "is_pathogenic": None,  # VUS
            "gene": gene,
        })

    # Shuffle to mix classes
    rng.shuffle(variants)

    # Write TSV
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    with open(output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["chr", "pos", "ref", "alt", "CLNSIG", "CLNREVSTAT", "GENEINFO"])
        for v in variants:
            writer.writerow([
                v["chr"], v["pos"], v["ref"], v["alt"],
                v["CLNSIG"], v["CLNREVSTAT"], v["GENEINFO"]
            ])

    print(f"Generated ClinVar data: {output}")
    print(f"  - Pathogenic: {n_pathogenic}")
    print(f"  - Benign: {n_benign}")
    print(f"  - VUS (excluded): {n_vus}")

    return variants


def generate_synthetic_dbnsfp(
    output_path: str,
    variants: list[dict],
    seed: int = 42,
) -> None:
    """
    Generate synthetic dbNSFP-like TSV data with prediction scores.

    Score characteristics:
    - CADD_phred: Good discriminator (AUC ~0.90), ~5% missing
    - REVEL_score: Excellent discriminator (AUC ~0.95), ~2% missing
    - MetaLR_score: Moderate discriminator (AUC ~0.75), ~10% missing
    - VEST4_score: Will be excluded due to ~40% missingness

    Args:
        output_path: Path to write the TSV file
        variants: List of variant dicts from generate_synthetic_clinvar
        seed: Random seed for reproducibility
    """
    rng = random.Random(seed + 1)  # Different seed for scores

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    with open(output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "chr", "pos", "ref", "alt",
            "CADD_phred", "REVEL_score", "MetaLR_score", "VEST4_score"
        ])

        for v in variants:
            is_path = v["is_pathogenic"]

            # For VUS, use random assignment for score generation
            if is_path is None:
                is_path = rng.choice([True, False])

            # CADD_phred: good discriminator, ~5% missing
            if rng.random() > 0.05:
                cadd = generate_score(is_path, 25, 5, 10, 4, 0, 50, rng)
                cadd_str = f"{cadd:.2f}"
            else:
                cadd_str = ""

            # REVEL_score: excellent discriminator, ~2% missing
            if rng.random() > 0.02:
                revel = generate_score(is_path, 0.75, 0.15, 0.20, 0.12, 0, 1, rng)
                revel_str = f"{revel:.4f}"
            else:
                revel_str = ""

            # MetaLR_score: moderate discriminator, ~10% missing
            if rng.random() > 0.10:
                metalr = generate_score(is_path, 0.65, 0.20, 0.35, 0.20, 0, 1, rng)
                metalr_str = f"{metalr:.4f}"
            else:
                metalr_str = ""

            # VEST4_score: ~40% missing (will be excluded)
            if rng.random() > 0.40:
                vest4 = generate_score(is_path, 0.60, 0.25, 0.40, 0.25, 0, 1, rng)
                vest4_str = f"{vest4:.4f}"
            else:
                vest4_str = ""

            writer.writerow([
                v["chr"], v["pos"], v["ref"], v["alt"],
                cadd_str, revel_str, metalr_str, vest4_str
            ])

    print(f"Generated dbNSFP data: {output}")


def generate_gene_list(output_path: str) -> None:
    """Generate a gene list file for testing --genes-file option."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    with open(output, "w") as f:
        for gene in GENE_INFO.keys():
            f.write(f"{gene}\n")

    print(f"Generated gene list: {output}")


def generate_variant_list(output_path: str, variants: list[dict]) -> None:
    """Generate a variant list file for testing --variants option."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    with open(output, "w") as f:
        f.write("# Test variants for PSROC\n")
        f.write("# Format: chr:pos:ref:alt\n")
        for v in variants:
            f.write(f"{v['chr']}:{v['pos']}:{v['ref']}:{v['alt']}\n")

    print(f"Generated variant list: {output}")


def build_hail_tables(
    clinvar_tsv: str,
    dbnsfp_tsv: str,
    output_dir: str,
    reference_genome: str = "GRCh38",
) -> tuple[str, str]:
    """
    Build Hail Tables from synthetic TSV files.

    This function requires Hail to be available. It creates properly
    structured Hail Tables that can be used with the PSROC pipeline.

    Args:
        clinvar_tsv: Path to synthetic ClinVar TSV
        dbnsfp_tsv: Path to synthetic dbNSFP TSV
        output_dir: Directory to write Hail Tables
        reference_genome: Reference genome (GRCh37 or GRCh38)

    Returns:
        Tuple of (clinvar_ht_path, dbnsfp_ht_path)
    """
    try:
        import hail as hl
        from hvantk.core.hail_context import init_hail
    except ImportError:
        raise ImportError("Hail is required to build Hail Tables. Install with: poetry install")

    init_hail()

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    # Build ClinVar table
    clinvar_ht_path = str(output / "synthetic_clinvar.ht")

    ht = hl.import_table(
        clinvar_tsv,
        types={
            "chr": hl.tstr,
            "pos": hl.tint32,
            "ref": hl.tstr,
            "alt": hl.tstr,
            "CLNSIG": hl.tstr,
            "CLNREVSTAT": hl.tstr,
            "GENEINFO": hl.tstr,
        },
    )

    # Create locus and alleles
    ht = ht.annotate(
        locus=hl.locus(ht.chr, ht.pos, reference_genome=reference_genome),
        alleles=hl.array([ht.ref, ht.alt]),
    )

    # Structure as ClinVar-like with info struct
    ht = ht.annotate(
        info=hl.struct(
            CLNSIG=hl.array([ht.CLNSIG]),
            CLNREVSTAT=hl.array([ht.CLNREVSTAT]),
            GENEINFO=hl.array([ht.GENEINFO]),
        )
    )

    ht = ht.select("locus", "alleles", "info")
    ht = ht.key_by("locus", "alleles")
    ht = ht.checkpoint(clinvar_ht_path, overwrite=True)

    print(f"Built ClinVar Hail Table: {clinvar_ht_path}")
    print(f"  - Variants: {ht.count()}")

    # Build dbNSFP table
    dbnsfp_ht_path = str(output / "synthetic_dbnsfp.ht")

    ht = hl.import_table(
        dbnsfp_tsv,
        types={
            "chr": hl.tstr,
            "pos": hl.tint32,
            "ref": hl.tstr,
            "alt": hl.tstr,
            "CADD_phred": hl.tstr,
            "REVEL_score": hl.tstr,
            "MetaLR_score": hl.tstr,
            "VEST4_score": hl.tstr,
        },
        missing="",
    )

    # Create locus and alleles
    ht = ht.annotate(
        locus=hl.locus(ht.chr, ht.pos, reference_genome=reference_genome),
        alleles=hl.array([ht.ref, ht.alt]),
    )

    # Convert score strings to floats
    ht = ht.annotate(
        CADD_phred=hl.if_else(ht.CADD_phred == "", hl.missing(hl.tfloat64), hl.float64(ht.CADD_phred)),
        REVEL_score=hl.if_else(ht.REVEL_score == "", hl.missing(hl.tfloat64), hl.float64(ht.REVEL_score)),
        MetaLR_score=hl.if_else(ht.MetaLR_score == "", hl.missing(hl.tfloat64), hl.float64(ht.MetaLR_score)),
        VEST4_score=hl.if_else(ht.VEST4_score == "", hl.missing(hl.tfloat64), hl.float64(ht.VEST4_score)),
    )

    ht = ht.select("locus", "alleles", "CADD_phred", "REVEL_score", "MetaLR_score", "VEST4_score")
    ht = ht.key_by("locus", "alleles")
    ht = ht.checkpoint(dbnsfp_ht_path, overwrite=True)

    print(f"Built dbNSFP Hail Table: {dbnsfp_ht_path}")
    print(f"  - Variants: {ht.count()}")

    return clinvar_ht_path, dbnsfp_ht_path


def main(output_dir: Optional[str] = None, seed: int = 42) -> None:
    """Generate all synthetic test data files."""
    if output_dir is None:
        output_dir = str(Path(__file__).parent)

    output = Path(output_dir)

    print("=" * 60)
    print("Generating PSROC Synthetic Test Data")
    print("=" * 60)

    # Generate ClinVar data
    clinvar_path = output / "synthetic_clinvar.tsv"
    variants = generate_synthetic_clinvar(str(clinvar_path), seed=seed)

    print()

    # Generate dbNSFP data
    dbnsfp_path = output / "synthetic_dbnsfp.tsv"
    generate_synthetic_dbnsfp(str(dbnsfp_path), variants, seed=seed)

    print()

    # Generate gene list
    genes_path = output / "test_genes.txt"
    generate_gene_list(str(genes_path))

    print()

    # Generate variant list
    variants_path = output / "test_variants.txt"
    generate_variant_list(str(variants_path), variants)

    print()
    print("=" * 60)
    print("Synthetic data generation complete!")
    print("=" * 60)
    print()
    print("To build Hail Tables (requires Hail):")
    print("  from generate_synthetic_data import build_hail_tables")
    print(f"  build_hail_tables('{clinvar_path}', '{dbnsfp_path}', '/tmp/psroc_test')")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate synthetic PSROC test data")
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: same directory as this script)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)",
    )
    args = parser.parse_args()

    main(output_dir=args.output_dir, seed=args.seed)