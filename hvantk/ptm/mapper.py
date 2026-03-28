"""Local protein-to-genomic coordinate mapper.

Maps PTM sites from protein coordinates (UniProt accession + residue position)
to genomic coordinates (chromosome + codon interval) using Ensembl GTF data.

Validated on 100 proteins in Phase 0.5 (99.1% mapping rate, 0.076s for 449 sites).
"""

import gzip
import logging
import re
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, NamedTuple

logger = logging.getLogger(__name__)


class CodonMapping(NamedTuple):
    """Result of mapping a protein residue to genomic coordinates."""

    chrom: str
    codon_start: int
    codon_end: int
    strand: str


class GTFData(NamedTuple):
    """Parsed GTF data for coordinate mapping."""

    cds_by_transcript: Dict[
        str, List[Tuple]
    ]  # ENST -> [(chrom, start, end, strand, phase)]
    mane_transcripts: set  # set of ENST IDs tagged as MANE Select
    gene_to_mane: Dict[str, str]  # gene_symbol -> ENST (MANE Select)
    transcript_to_gene: Dict[str, str]  # ENST -> gene_symbol


def parse_ensembl_gtf(gtf_path: str) -> GTFData:
    """Parse Ensembl GTF in a single pass to extract CDS exons, MANE Select tags, and gene names.

    Parameters
    ----------
    gtf_path : str
        Path to the Ensembl GTF file (gzipped or plain text).

    Returns
    -------
    GTFData
        Parsed data with CDS exon intervals, MANE Select tags, and gene mappings.
    """
    cds_by_transcript = defaultdict(list)
    mane_transcripts = set()
    gene_to_mane = {}
    transcript_to_gene = {}

    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            if fields[2] == "transcript":
                m_tid = re.search(r'transcript_id "([^"]+)"', fields[8])
                m_gene = re.search(r'gene_name "([^"]+)"', fields[8])
                if m_tid and m_gene:
                    enst = m_tid.group(1).split(".")[0]
                    gene = m_gene.group(1)
                    transcript_to_gene[enst] = gene
                    if 'tag "MANE_Select"' in fields[8]:
                        mane_transcripts.add(enst)
                        gene_to_mane[gene] = enst

            if fields[2] != "CDS":
                continue

            chrom = fields[0]
            start = int(fields[3])  # 1-based inclusive
            end = int(fields[4])  # 1-based inclusive
            strand = fields[6]
            phase = int(fields[7])

            m = re.search(r'transcript_id "([^"]+)"', fields[8])
            if not m:
                continue
            enst = m.group(1).split(".")[0]
            cds_by_transcript[enst].append((chrom, start, end, strand, phase))

    # Sort exons by genomic position within each transcript
    for enst in cds_by_transcript:
        cds_by_transcript[enst].sort(key=lambda x: x[1])

    logger.info(
        f"Parsed GTF: {sum(len(v) for v in cds_by_transcript.values()):,} CDS features, "
        f"{len(cds_by_transcript):,} transcripts, "
        f"{len(mane_transcripts):,} MANE Select"
    )

    return GTFData(
        cds_by_transcript=dict(cds_by_transcript),
        mane_transcripts=mane_transcripts,
        gene_to_mane=gene_to_mane,
        transcript_to_gene=transcript_to_gene,
    )


def map_residue_to_genomic(
    enst_id: str,
    residue_pos: int,
    cds_lookup: Dict[str, List[Tuple]],
) -> Optional[CodonMapping]:
    """Map a single protein residue to its genomic codon coordinates.

    Parameters
    ----------
    enst_id : str
        Ensembl transcript ID (without version suffix).
    residue_pos : int
        1-based protein residue position.
    cds_lookup : dict
        Mapping of ENST -> sorted list of (chrom, start, end, strand, phase).

    Returns
    -------
    CodonMapping or None
        Genomic coordinates of the codon, or None if mapping fails.
    """
    exons = cds_lookup.get(enst_id)
    if not exons:
        return None

    chrom = exons[0][0]
    strand = exons[0][3]

    # Order exons in CDS reading direction (5'->3' of mRNA)
    if strand == "+":
        ordered = sorted(exons, key=lambda x: x[1])
    else:
        ordered = sorted(exons, key=lambda x: x[1], reverse=True)

    # Build flat genomic position list in CDS order
    positions = []
    for _, start, end, s, _ in ordered:
        if s == "+":
            positions.extend(range(start, end + 1))
        else:
            positions.extend(range(end, start - 1, -1))

    # Get 3 positions for the codon
    idx = (residue_pos - 1) * 3
    if idx + 3 > len(positions):
        return None  # out-of-bounds: protein longer than CDS

    codon_pos = positions[idx : idx + 3]
    return CodonMapping(
        chrom=chrom,
        codon_start=min(codon_pos),
        codon_end=max(codon_pos),
        strand=strand,
    )


def map_protein_sites(
    enst_id: str,
    residue_positions: List[int],
    cds_lookup: Dict[str, List[Tuple]],
) -> List[Tuple[int, Optional[CodonMapping]]]:
    """Batch-map multiple residue positions for one protein.

    Builds the position list once and reuses it for all positions.

    Parameters
    ----------
    enst_id : str
        Ensembl transcript ID (without version suffix).
    residue_positions : list of int
        1-based protein residue positions to map.
    cds_lookup : dict
        Mapping of ENST -> sorted list of (chrom, start, end, strand, phase).

    Returns
    -------
    list of (int, CodonMapping or None)
        Tuples of (residue_pos, mapping_or_None).
    """
    exons = cds_lookup.get(enst_id)
    if not exons:
        return [(p, None) for p in residue_positions]

    chrom = exons[0][0]
    strand = exons[0][3]

    if strand == "+":
        ordered = sorted(exons, key=lambda x: x[1])
    else:
        ordered = sorted(exons, key=lambda x: x[1], reverse=True)

    positions = []
    for _, start, end, s, _ in ordered:
        if s == "+":
            positions.extend(range(start, end + 1))
        else:
            positions.extend(range(end, start - 1, -1))

    results = []
    for p in residue_positions:
        idx = (p - 1) * 3
        if idx + 3 > len(positions):
            results.append((p, None))
        else:
            codon_pos = positions[idx : idx + 3]
            results.append(
                (
                    p,
                    CodonMapping(
                        chrom=chrom,
                        codon_start=min(codon_pos),
                        codon_end=max(codon_pos),
                        strand=strand,
                    ),
                )
            )
    return results


def resolve_transcript(
    ensembl_xrefs: List[Dict],
    gene_symbol: str,
    gtf_data: GTFData,
) -> Tuple[Optional[str], str]:
    """Resolve UniProt Ensembl cross-refs to a single transcript ID.

    Uses the 3-strategy cascade validated in Phase 0.5:
    1. MANE Select xref (97% of proteins in Phase 0.5)
    2. Any Ensembl xref with CDS data (3%)
    3. Gene name -> MANE Select fallback (0% needed, safety net)

    Parameters
    ----------
    ensembl_xrefs : list of dict
        UniProt Ensembl cross-references (each has 'id' key).
    gene_symbol : str
        Gene symbol for fallback resolution.
    gtf_data : GTFData
        Parsed GTF data.

    Returns
    -------
    (str or None, str)
        Tuple of (ENST ID or None, resolution method name).
    """
    # Strategy 1: Prefer MANE Select
    for xref in ensembl_xrefs:
        candidate = xref.get("id", "").split(".")[0]
        if (
            candidate in gtf_data.mane_transcripts
            and candidate in gtf_data.cds_by_transcript
        ):
            return candidate, "xref_mane"

    # Strategy 2: Any xref with CDS data
    for xref in ensembl_xrefs:
        candidate = xref.get("id", "").split(".")[0]
        if candidate in gtf_data.cds_by_transcript:
            return candidate, "xref_any"

    # Strategy 3: Gene name -> MANE Select fallback
    if gene_symbol in gtf_data.gene_to_mane:
        candidate = gtf_data.gene_to_mane[gene_symbol]
        if candidate in gtf_data.cds_by_transcript:
            return candidate, "gene_mane"

    return None, "unresolved"
