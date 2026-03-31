"""
PTM Build Pipeline — orchestrates PTM data acquisition, coordinate mapping, and table building.

This module exposes the build workflow as a Python API so it can be used
programmatically (notebooks, scripts) or from the CLI.

Example:
    >>> from hvantk.ptm.pipeline import PTMBuildConfig, ptm_build_pipeline
    >>> config = PTMBuildConfig(
    ...     output_dir="data/ptm/",
    ...     output_ht="data/ptm/ptm_sites.ht",
    ... )
    >>> result = ptm_build_pipeline(config)
    >>> print(result.n_mapped, result.mapped_tsv_path)
"""

import csv
import logging
import os
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from hvantk.ptm.constants import (
    ENSEMBL_GTF_URL,
    ENSEMBL_GTF_FILENAME,
    PTM_TYPE_CATEGORIES,
    PTM_OUTPUT_COLUMNS,
)
from hvantk.ptm.mapper import (
    GTFData,
    parse_ensembl_gtf,
    map_protein_sites,
    resolve_transcript,
)

logger = logging.getLogger(__name__)


@dataclass
class PTMBuildConfig:
    """Configuration for the PTM build pipeline.

    Attributes:
        output_dir: Directory for intermediate files (GTF, UniProt TSV, mapped TSV).
        output_ht: Path to write the final PTM sites Hail Table.
        gtf_path: Path to a pre-downloaded Ensembl GTF (None = auto-download).
        ptm_tsv: Path to a pre-downloaded UniProt PTM TSV (None = auto-download).
        flanking_codons: Number of flanking codons for proximal window.
        reference_genome: Reference genome for Hail Table.
        overwrite: Whether to overwrite existing outputs.
    """

    output_dir: str = ""
    output_ht: str = ""
    gtf_path: Optional[str] = None
    ptm_tsv: Optional[str] = None
    peptideatlas_tsv: Optional[str] = None
    cptac_tsv: Optional[str] = None
    flanking_codons: int = 5
    reference_genome: str = "GRCh38"
    overwrite: bool = False

    def validate(self) -> List[str]:
        """Validate configuration and return list of errors."""
        errors = []
        if not self.output_dir:
            errors.append("output_dir is required")
        if not self.output_ht:
            errors.append("output_ht is required")
        if self.gtf_path and not os.path.exists(self.gtf_path):
            errors.append(f"GTF file not found: {self.gtf_path}")
        if self.ptm_tsv and not os.path.exists(self.ptm_tsv):
            errors.append(f"PTM TSV file not found: {self.ptm_tsv}")
        if self.peptideatlas_tsv and not os.path.exists(self.peptideatlas_tsv):
            errors.append(f"PeptideAtlas TSV file not found: {self.peptideatlas_tsv}")
        if self.cptac_tsv and not os.path.exists(self.cptac_tsv):
            errors.append(f"CPTAC TSV file not found: {self.cptac_tsv}")
        return errors


@dataclass
class PTMBuildResult:
    """Result from a PTM build pipeline run.

    Attributes:
        n_total: Total PTM sites processed.
        n_mapped: Number of sites successfully mapped to genomic coordinates.
        n_failed: Number of sites that failed mapping.
        resolution_counts: Count per transcript resolution method.
        mapped_tsv_path: Path to the mapped TSV file.
        output_ht: Path to the built Hail Table.
        gtf_stats: GTF parsing statistics.
    """

    n_total: int = 0
    n_mapped: int = 0
    n_failed: int = 0
    resolution_counts: Dict[str, int] = field(default_factory=dict)
    mapped_tsv_path: str = ""
    output_ht: str = ""
    gtf_stats: Dict[str, int] = field(default_factory=dict)


def download_ensembl_gtf(output_dir: str, overwrite: bool = False) -> str:
    """Download the Ensembl GTF file if not already cached.

    Parameters
    ----------
    output_dir : str
        Directory to save the GTF file.
    overwrite : bool
        If True, re-download even if the file exists.

    Returns
    -------
    str
        Path to the downloaded GTF file.
    """
    os.makedirs(output_dir, exist_ok=True)
    gtf_path = os.path.join(output_dir, ENSEMBL_GTF_FILENAME)

    if os.path.exists(gtf_path) and not overwrite:
        logger.info(f"Using cached GTF: {gtf_path}")
        return gtf_path

    logger.info(f"Downloading Ensembl GTF to {gtf_path}...")
    from urllib.request import urlretrieve

    urlretrieve(ENSEMBL_GTF_URL, gtf_path)
    logger.info(f"Downloaded: {gtf_path}")
    return gtf_path


def download_uniprot_ptm(output_dir: str, overwrite: bool = False) -> str:
    """Download UniProt PTM data via the REST API.

    Parameters
    ----------
    output_dir : str
        Directory to save the TSV file.
    overwrite : bool
        If True, re-download even if a file exists.

    Returns
    -------
    str
        Path to the downloaded TSV file.
    """
    from hvantk.datasets.uniprot_ptm_datasets import UniProtPTMDataset

    dataset = UniProtPTMDataset.from_latest()
    return dataset.download(output_dir, overwrite=overwrite)


def map_ptm_sites(
    ptm_tsv: str,
    gtf_data: GTFData,
    output_path: str,
) -> PTMBuildResult:
    """Map PTM sites from protein coordinates to genomic coordinates.

    Reads the UniProt PTM TSV, resolves transcripts, maps each site to
    genomic codon coordinates, and writes the result to a new TSV.

    Parameters
    ----------
    ptm_tsv : str
        Path to the UniProt PTM TSV (from download_uniprot_ptm).
    gtf_data : GTFData
        Parsed Ensembl GTF data.
    output_path : str
        Path to write the mapped TSV.

    Returns
    -------
    PTMBuildResult
        Mapping statistics.
    """
    resolution_counts: Dict[str, int] = defaultdict(int)
    n_mapped = 0
    n_failed = 0
    n_total = 0

    with open(ptm_tsv) as fin, open(output_path, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(
            fout, fieldnames=PTM_OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()

        current_acc = None
        current_records: List[dict] = []

        def flush_protein(records: List[dict]) -> None:
            nonlocal n_mapped, n_failed, n_total
            if not records:
                return

            xref_str = records[0].get("ensembl_xrefs", "")
            xrefs = [{"id": xid.strip()} for xid in xref_str.split(";") if xid.strip()]
            gene = records[0].get("gene_symbol", "")

            enst, method = resolve_transcript(xrefs, gene, gtf_data)
            resolution_counts[method] += 1

            if enst is None:
                n_failed += len(records)
                n_total += len(records)
                return

            positions = [int(r["position"]) for r in records]
            mappings = map_protein_sites(enst, positions, gtf_data.cds_by_transcript)
            pos_to_mapping = dict(mappings)

            for rec in records:
                n_total += 1
                pos = int(rec["position"])
                m = pos_to_mapping.get(pos)
                if m is None:
                    n_failed += 1
                    continue

                desc = rec.get("description", "")
                ptm_category = "other"
                for prefix, cat in PTM_TYPE_CATEGORIES.items():
                    if desc.startswith(prefix):
                        ptm_category = cat
                        break

                writer.writerow(
                    {
                        "chrom": m.chrom,
                        "codon_start": m.codon_start,
                        "codon_end": m.codon_end,
                        "strand": m.strand,
                        "uniprot_id": rec.get("accession", ""),
                        "gene_symbol": gene,
                        "residue_pos": pos,
                        "amino_acid": rec.get("amino_acid", ""),
                        "ptm_type": desc,
                        "ptm_category": ptm_category,
                        "source_db": rec.get("source_db", "UniProt"),
                        "evidence_type": rec.get("evidence_type", "curated"),
                        "n_observations": rec.get("n_observations", "0"),
                    }
                )
                n_mapped += 1

        for row in reader:
            acc = row.get("accession", "")
            if acc != current_acc:
                flush_protein(current_records)
                current_acc = acc
                current_records = []
            current_records.append(row)

        flush_protein(current_records)

    logger.info(
        f"Mapping complete: {n_mapped}/{n_total} mapped "
        f"({100 * n_mapped / max(n_total, 1):.1f}%), {n_failed} failed"
    )
    logger.info(f"Resolution: {dict(resolution_counts)}")

    return PTMBuildResult(
        n_total=n_total,
        n_mapped=n_mapped,
        n_failed=n_failed,
        resolution_counts=dict(resolution_counts),
        mapped_tsv_path=output_path,
    )


def ptm_build_pipeline(config: PTMBuildConfig) -> PTMBuildResult:
    """Run the full PTM build pipeline.

    Steps:
        1. Download Ensembl GTF (if not provided)
        2. Download UniProt PTM data (if not provided)
        3. Parse GTF and map PTM sites to genomic coordinates
        4. Build Hail Table from mapped coordinates

    Parameters
    ----------
    config : PTMBuildConfig
        Pipeline configuration.

    Returns
    -------
    PTMBuildResult
        Mapping statistics and output paths.

    Raises
    ------
    ValueError
        If configuration validation fails.
    """
    errors = config.validate()
    if errors:
        raise ValueError(f"Invalid config: {'; '.join(errors)}")

    os.makedirs(config.output_dir, exist_ok=True)

    # Step 1: Ensure GTF is available
    gtf_path = config.gtf_path or download_ensembl_gtf(
        config.output_dir, config.overwrite
    )

    # Step 2: Parse GTF
    logger.info("Parsing Ensembl GTF...")
    gtf_data = parse_ensembl_gtf(gtf_path)

    # Step 3: Ensure UniProt PTM TSV is available
    ptm_tsv = config.ptm_tsv
    if ptm_tsv is None:
        ptm_tsv = download_uniprot_ptm(config.output_dir, config.overwrite)

    # Step 4: Map PTM sites to genomic coordinates
    mapped_path = os.path.join(config.output_dir, "ptm_sites_mapped.tsv")
    result = map_ptm_sites(ptm_tsv, gtf_data, mapped_path)

    # Step 4b: Map PeptideAtlas sites (if provided)
    if config.peptideatlas_tsv:
        logger.info("Mapping PeptideAtlas phospho sites...")
        pa_mapped_path = os.path.join(config.output_dir, "peptideatlas_sites_mapped.tsv")
        pa_result = map_ptm_sites(config.peptideatlas_tsv, gtf_data, pa_mapped_path)

        # Concatenate mapped TSVs
        combined_path = os.path.join(config.output_dir, "ptm_sites_combined.tsv")
        with open(combined_path, "w", newline="") as fout:
            writer = csv.DictWriter(
                fout, fieldnames=PTM_OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for src_path in [mapped_path, pa_mapped_path]:
                with open(src_path) as fin:
                    reader = csv.DictReader(fin, delimiter="\t")
                    for row in reader:
                        writer.writerow(row)

        mapped_path = combined_path
        result.n_total += pa_result.n_total
        result.n_mapped += pa_result.n_mapped
        result.n_failed += pa_result.n_failed
        for method, count in pa_result.resolution_counts.items():
            result.resolution_counts[method] = (
                result.resolution_counts.get(method, 0) + count
            )

        logger.info(
            f"Combined: {result.n_mapped} total mapped sites "
            f"(UniProt + PeptideAtlas)"
        )

    # Step 4c: Map CPTAC sites (if provided)
    if config.cptac_tsv:
        logger.info("Mapping CPTAC phospho sites...")
        cptac_mapped_path = os.path.join(config.output_dir, "cptac_sites_mapped.tsv")
        cptac_result = map_ptm_sites(config.cptac_tsv, gtf_data, cptac_mapped_path)

        # Concatenate with existing mapped TSV
        combined_path = os.path.join(config.output_dir, "ptm_sites_combined.tsv")
        with open(combined_path, "w", newline="") as fout:
            writer = csv.DictWriter(
                fout, fieldnames=PTM_OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for src_path in [mapped_path, cptac_mapped_path]:
                with open(src_path) as fin:
                    reader = csv.DictReader(fin, delimiter="\t")
                    for row in reader:
                        writer.writerow(row)

        mapped_path = combined_path
        result.n_total += cptac_result.n_total
        result.n_mapped += cptac_result.n_mapped
        result.n_failed += cptac_result.n_failed
        for method, count in cptac_result.resolution_counts.items():
            result.resolution_counts[method] = (
                result.resolution_counts.get(method, 0) + count
            )

        logger.info(
            f"Combined: {result.n_mapped} total mapped sites "
            f"(including CPTAC)"
        )

    # Step 5: Build Hail Table
    logger.info(f"Building Hail Table at {config.output_ht}...")
    from hvantk.tables.table_builders import create_ptm_sites_tb

    create_ptm_sites_tb(
        input_path=mapped_path,
        output_path=config.output_ht,
        reference_genome=config.reference_genome,
        flanking_codons=config.flanking_codons,
        overwrite=config.overwrite,
    )

    result.output_ht = config.output_ht
    result.gtf_stats = {
        "transcripts": len(gtf_data.cds_by_transcript),
        "mane_select": len(gtf_data.mane_transcripts),
        "genes_with_mane": len(gtf_data.gene_to_mane),
    }

    logger.info(f"PTM build pipeline complete: {result.n_mapped} sites mapped")
    return result
