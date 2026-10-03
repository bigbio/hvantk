"""
PTM Build Pipeline — pure coordinate-mapping core for PTM data acquisition.

This module exposes the algorithm-only (no skills) build core as a Python API.
Skill-driven steps (UniProt download, Hail Table build) live in
:mod:`hvantk.tools.ptm.pipeline`.

Example:
    >>> from hvantk.algorithms.ptm.pipeline import PTMBuildConfig, ptm_build_pipeline_core
    >>> config = PTMBuildConfig(
    ...     output_dir="data/ptm/",
    ...     ptm_tsv="data/ptm/uniprot-ptm-human-2026-10-01.tsv",
    ... )
    >>> result = ptm_build_pipeline_core(config)
    >>> print(result.source_counts, result.n_mapped, result.mapped_tsv_path)
"""

import csv
import gzip
import logging
import os
import shutil
import urllib.parse

import requests
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from hvantk.core.utils.bgzf import BgzfWriter
from hvantk.algorithms.ptm.constants import (
    ENSEMBL_GTF_URL,
    ENSEMBL_GTF_FILENAME,
    PTM_TYPE_CATEGORIES,
    PTM_OUTPUT_COLUMNS,
)
from hvantk.algorithms.ptm.mapper import (
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
        output_ht: Path for the PTM sites Hail Table. Only
            :func:`hvantk.tools.ptm.pipeline.ptm_build_pipeline` writes it, and it
            requires it; :func:`ptm_build_pipeline_core` ignores it.
        gtf_path: Path to a pre-downloaded Ensembl GTF (None = auto-download).
        ptm_tsv: Path to a pre-downloaded UniProt PTM TSV. Required by
            :func:`ptm_build_pipeline_core`; when None,
            :func:`hvantk.tools.ptm.pipeline.ptm_build_pipeline` downloads it.
        peptideatlas_tsv: PeptideAtlas phospho TSV (``hvantk download
            peptideatlas-phospho``); its sites are added when set.
        cptac_tsv: CPTAC phospho TSV (``hvantk download cptac-phospho``); its sites
            are added when set.
        flanking_codons: Number of flanking codons for the Hail Table's proximal
            window (``flanking_interval``); :func:`ptm_build_pipeline_core` ignores it.
        reference_genome: Reference genome for the Hail Table;
            :func:`ptm_build_pipeline_core` ignores it.
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
        """Return the errors that stop every caller; empty when there are none.

        :func:`ptm_build_pipeline_core` also requires ``ptm_tsv``, and
        :func:`hvantk.tools.ptm.pipeline.ptm_build_pipeline` also requires
        ``output_ht``; each checks its own field before doing any work.
        """
        errors = []
        if not self.output_dir:
            errors.append("output_dir is required")
        if self.flanking_codons < 0:
            errors.append(f"flanking_codons must be >= 0, got {self.flanking_codons}")
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
        mapped_tsv_path: Path to the mapped TSV file (the combined TSV when more
            than one source was mapped).
        output_ht: Path to the built Hail Table; empty when none was built (the
            core never builds one, and the tools pipeline skips it when no site maps).
        gtf_stats: GTF parsing statistics.
        source_counts: Sites mapped per source, in mapping order: "UniProt", then
            "PeptideAtlas" and/or "CPTAC" when their TSVs were set. Each source is
            named after the config field that held its file, not read from the
            rows' ``source_db``, and can be 0 (the core then logs a warning).
    """

    n_total: int = 0
    n_mapped: int = 0
    n_failed: int = 0
    resolution_counts: Dict[str, int] = field(default_factory=dict)
    mapped_tsv_path: str = ""
    output_ht: str = ""
    gtf_stats: Dict[str, int] = field(default_factory=dict)
    source_counts: Dict[str, int] = field(default_factory=dict)


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
        logger.info("Using cached GTF: %s", gtf_path)
        return gtf_path

    parsed = urllib.parse.urlparse(ENSEMBL_GTF_URL)
    if parsed.scheme != "https":
        raise ValueError(
            "Invalid Ensembl GTF URL scheme (expected https): %s" % ENSEMBL_GTF_URL
        )
    host = parsed.hostname or ""
    if host != "ftp.ensembl.org" and not host.endswith(".ensembl.org"):
        raise ValueError(
            "Invalid Ensembl GTF URL host "
            "(expected trusted Ensembl host, e.g. ftp.ensembl.org or *.ensembl.org): %s"
            % ENSEMBL_GTF_URL
        )

    logger.info("Downloading Ensembl GTF to %s...", gtf_path)
    with requests.get(ENSEMBL_GTF_URL, stream=True, timeout=600) as resp:
        resp.raise_for_status()
        with open(gtf_path, "wb") as fout:
            shutil.copyfileobj(resp.raw, fout)
    logger.info("Downloaded: %s", gtf_path)
    return gtf_path


def map_ptm_sites(
    ptm_tsv: str,
    gtf_data: GTFData,
    output_path: str,
    transcript_cache: Optional[Dict] = None,
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
    transcript_cache : dict, optional
        Shared TranscriptCDS cache for cross-call reuse. If None, a local
        cache is created for this call only.

    Returns
    -------
    PTMBuildResult
        Mapping statistics.
    """
    if transcript_cache is None:
        transcript_cache = {}

    resolution_counts: Dict[str, int] = defaultdict(int)
    n_mapped = 0
    n_failed = 0
    n_total = 0

    # Write BGZF so Hail can import in parallel across Spark partitions.
    if not output_path.endswith((".bgz", ".gz")):
        output_path += ".bgz"

    with open(ptm_tsv, encoding="utf-8") as fin, BgzfWriter(output_path) as fout:
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
            mappings = map_protein_sites(
                enst, positions, gtf_data.cds_by_transcript, cache=transcript_cache
            )
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
                        "tissue_type": rec.get("tissue_type", ""),
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
        "Mapping complete: %d/%d mapped (%.1f%%), %d failed",
        n_mapped,
        n_total,
        100 * n_mapped / max(n_total, 1),
        n_failed,
    )
    logger.info("Resolution: %s", dict(resolution_counts))

    return PTMBuildResult(
        n_total=n_total,
        n_mapped=n_mapped,
        n_failed=n_failed,
        resolution_counts=dict(resolution_counts),
        mapped_tsv_path=output_path,
    )


def _concat_bgz_tsvs(src_paths: List[str], dest_path: str) -> None:
    """Concatenate multiple BGZF TSV files into a single BGZF TSV.

    Reads each source with ``gzip.open`` (which handles BGZF transparently),
    writes a single header from the first file, and streams data rows into a
    new BGZF file via :class:`BgzfWriter`.
    """
    with BgzfWriter(dest_path) as fout:
        writer = None
        for src_path in src_paths:
            with gzip.open(src_path, "rt", encoding="utf-8") as fin:
                reader = csv.DictReader(fin, delimiter="\t")
                if writer is None:
                    writer = csv.DictWriter(
                        fout,
                        fieldnames=reader.fieldnames,
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                for row in reader:
                    writer.writerow(row)


def ptm_build_pipeline_core(config: PTMBuildConfig) -> PTMBuildResult:
    """Run the pure coordinate-mapping core of the PTM build pipeline.

    Steps:
        1. Download Ensembl GTF (if not provided)
        2. Parse the GTF
        3. Map the UniProt sites to genomic coordinates
        3b. Map PeptideAtlas / CPTAC sites when their TSVs are set, and
            concatenate every source into ``ptm_sites_combined.tsv.bgz``

    The UniProt download and Hail Table build steps require
    :mod:`hvantk.skills` and are handled by the workflow layer in
    :mod:`hvantk.tools.ptm.pipeline`.

    Parameters
    ----------
    config : PTMBuildConfig
        Pipeline configuration.  ``config.ptm_tsv`` must be set to a
        pre-downloaded UniProt TSV path; callers in the tools layer are
        responsible for downloading it first when it is None.
        ``config.output_ht``, ``config.flanking_codons`` and
        ``config.reference_genome`` are not used here.

    Returns
    -------
    PTMBuildResult
        Mapping statistics, output paths (``mapped_tsv_path`` is set) and the
        sites mapped per source (``source_counts``).

    Raises
    ------
    ValueError
        If configuration validation fails or ``config.ptm_tsv`` is not set.
    """
    errors = config.validate()
    if errors:
        raise ValueError(f"Invalid config: {'; '.join(errors)}")
    # Checked before the GTF download and parse, so a missing UniProt TSV fails at once.
    if not config.ptm_tsv:
        raise ValueError(
            "config.ptm_tsv must be set before calling ptm_build_pipeline_core. "
            "Download it with hvantk.tools.ptm.pipeline.download_uniprot_ptm "
            "(`hvantk download uniprot-ptm`), or use "
            "hvantk.tools.ptm.pipeline.ptm_build_pipeline, which downloads it and "
            "also builds the Hail Table."
        )

    os.makedirs(config.output_dir, exist_ok=True)

    # Step 1: Ensure GTF is available
    gtf_path = config.gtf_path or download_ensembl_gtf(
        config.output_dir, config.overwrite
    )

    # Step 2: Parse GTF
    logger.info("Parsing Ensembl GTF...")
    gtf_data = parse_ensembl_gtf(gtf_path)

    # Shared transcript cache across all mapping calls — avoids rebuilding
    # TranscriptCDS objects when the same ENST appears in multiple sources.
    transcript_cache: Dict = {}

    # Step 3: Map the UniProt sites to genomic coordinates
    mapped_path = os.path.join(config.output_dir, "ptm_sites_mapped.tsv.bgz")
    result = map_ptm_sites(config.ptm_tsv, gtf_data, mapped_path, transcript_cache)
    # Taken now, before the extra sources are added into result.n_mapped.
    source_counts = {"UniProt": result.n_mapped}
    source_inputs = {"UniProt": config.ptm_tsv}

    # Step 3b: Map additional PTM sources (PeptideAtlas, CPTAC) and
    # merge everything into a single canonical combined TSV.
    extra_sources: List[tuple] = []
    if config.peptideatlas_tsv:
        extra_sources.append(
            (
                "PeptideAtlas",
                config.peptideatlas_tsv,
                "peptideatlas_sites_mapped.tsv.bgz",
            )
        )
    if config.cptac_tsv:
        extra_sources.append(("CPTAC", config.cptac_tsv, "cptac_sites_mapped.tsv.bgz"))

    source_paths = [mapped_path]
    for source_name, source_tsv, source_filename in extra_sources:
        logger.info("Mapping %s phospho sites...", source_name)
        source_mapped_path = os.path.join(config.output_dir, source_filename)
        source_result = map_ptm_sites(
            source_tsv, gtf_data, source_mapped_path, transcript_cache
        )
        source_paths.append(source_mapped_path)
        source_counts[source_name] = source_result.n_mapped
        source_inputs[source_name] = source_tsv
        result.n_total += source_result.n_total
        result.n_mapped += source_result.n_mapped
        result.n_failed += source_result.n_failed
        for method, count in source_result.resolution_counts.items():
            result.resolution_counts[method] = (
                result.resolution_counts.get(method, 0) + count
            )

    result.source_counts = source_counts
    for name, n_mapped in source_counts.items():
        if n_mapped == 0:
            # Usually the wrong file: say which one, since the totals hide it.
            logger.warning(
                "%s: none of the PTM sites in %s mapped to the genome, so it adds "
                "nothing to the output",
                name,
                source_inputs[name],
            )

    if len(source_paths) > 1:
        combined_path = os.path.join(config.output_dir, "ptm_sites_combined.tsv.bgz")
        _concat_bgz_tsvs(source_paths, combined_path)
        mapped_path = combined_path
        result.mapped_tsv_path = combined_path
        logger.info(
            "Combined: %d total mapped sites across %d sources",
            result.n_mapped,
            len(source_paths),
        )

    result.gtf_stats = {
        "transcripts": len(gtf_data.cds_by_transcript),
        "mane_select": len(gtf_data.mane_transcripts),
        "genes_with_mane": len(gtf_data.gene_to_mane),
    }

    logger.info("PTM core mapping complete: %d sites mapped", result.n_mapped)
    return result
