"""Local protein-to-genomic coordinate mapper.

Maps PTM sites from protein coordinates (UniProt accession + residue position)
to genomic coordinates (chromosome + codon interval) using Ensembl GTF data.

Validated on 100 proteins in Phase 0.5 (99.1% mapping rate, 0.076s for 449 sites).
"""

import bisect
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


class TranscriptCDS:
    """Pre-computed CDS structure for O(log n) residue-to-genomic mapping.

    Instead of materializing a list of every CDS genomic position (O(CDS length)),
    this stores cumulative exon offsets and uses binary search to resolve a CDS
    offset to a genomic coordinate in O(log exons) time.
    """

    __slots__ = ("chrom", "strand", "_cum_starts", "_exon_genomic", "_total_length")

    def __init__(self, exons: List[Tuple]) -> None:
        self.chrom: str = exons[0][0]
        self.strand: str = exons[0][3]

        if self.strand == "+":
            ordered = sorted(exons, key=lambda x: x[1])
        else:
            ordered = sorted(exons, key=lambda x: x[1], reverse=True)

        cum_starts: List[int] = []
        exon_genomic: List[Tuple[int, int]] = []
        cum = 0
        for _, start, end, s, _ in ordered:
            cum_starts.append(cum)
            length = end - start + 1
            if s == "+":
                exon_genomic.append((start, 1))
            else:
                exon_genomic.append((end, -1))
            cum += length

        self._cum_starts = cum_starts
        self._exon_genomic = exon_genomic
        self._total_length = cum

    def _cds_offset_to_genomic(self, cds_offset: int) -> int:
        """Map a 0-based CDS offset to a 1-based genomic position."""
        exon_idx = bisect.bisect_right(self._cum_starts, cds_offset) - 1
        offset_in_exon = cds_offset - self._cum_starts[exon_idx]
        base, direction = self._exon_genomic[exon_idx]
        return base + direction * offset_in_exon

    def map_residue(self, residue_pos: int) -> Optional[CodonMapping]:
        """Map a 1-based residue position to genomic codon coordinates."""
        idx = (residue_pos - 1) * 3
        if idx + 3 > self._total_length:
            return None
        p0 = self._cds_offset_to_genomic(idx)
        p1 = self._cds_offset_to_genomic(idx + 1)
        p2 = self._cds_offset_to_genomic(idx + 2)
        return CodonMapping(
            chrom=self.chrom,
            codon_start=min(p0, p1, p2),
            codon_end=max(p0, p1, p2),
            strand=self.strand,
        )


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

    # Pre-compile regexes — avoids recompilation on every line (~3.5M lines)
    _re_tid = re.compile(r'transcript_id "([^"]+)"')
    _re_gene = re.compile(r'gene_name "([^"]+)"')

    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue

            # Quick scan for feature type before full split — skip the ~80%
            # of lines that are gene/exon/UTR/start_codon/stop_codon/etc.
            # GTF columns are tab-separated; feature type is in column 3.
            tab1 = line.index("\t")
            tab2 = line.index("\t", tab1 + 1)
            tab3 = line.index("\t", tab2 + 1)
            feature = line[tab2 + 1 : tab3]

            if feature != "transcript" and feature != "CDS":
                continue

            # Only split for transcript + CDS lines (~5% of total)
            fields = line.split("\t", 9)
            if len(fields) < 9:
                continue
            attrs = fields[8]

            if feature == "transcript":
                m_tid = _re_tid.search(attrs)
                m_gene = _re_gene.search(attrs)
                if m_tid and m_gene:
                    enst = m_tid.group(1).split(".")[0]
                    gene = m_gene.group(1)
                    transcript_to_gene[enst] = gene
                    if 'tag "MANE_Select"' in attrs:
                        mane_transcripts.add(enst)
                        gene_to_mane[gene] = enst
                continue

            # CDS line
            chrom = fields[0]
            start = int(fields[3])  # 1-based inclusive
            end = int(fields[4])  # 1-based inclusive
            strand = fields[6]
            phase = int(fields[7])

            m = _re_tid.search(attrs)
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


def _get_transcript_cds(
    enst_id: str,
    cds_lookup: Dict[str, List[Tuple]],
    cache: Optional[Dict[str, TranscriptCDS]] = None,
) -> Optional[TranscriptCDS]:
    """Get or create a cached TranscriptCDS for the given transcript."""
    if cache is not None and enst_id in cache:
        return cache[enst_id]
    exons = cds_lookup.get(enst_id)
    if not exons:
        return None
    tcds = TranscriptCDS(exons)
    if cache is not None:
        cache[enst_id] = tcds
    return tcds


def map_residue_to_genomic(
    enst_id: str,
    residue_pos: int,
    cds_lookup: Dict[str, List[Tuple]],
    cache: Optional[Dict[str, TranscriptCDS]] = None,
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
    cache : dict, optional
        Shared TranscriptCDS cache for cross-call reuse.

    Returns
    -------
    CodonMapping or None
        Genomic coordinates of the codon, or None if mapping fails.
    """
    tcds = _get_transcript_cds(enst_id, cds_lookup, cache)
    if tcds is None:
        return None
    return tcds.map_residue(residue_pos)


def map_protein_sites(
    enst_id: str,
    residue_positions: List[int],
    cds_lookup: Dict[str, List[Tuple]],
    cache: Optional[Dict[str, TranscriptCDS]] = None,
) -> List[Tuple[int, Optional[CodonMapping]]]:
    """Batch-map multiple residue positions for one protein.

    Uses a pre-computed TranscriptCDS with binary search (O(log exons) per
    residue) instead of materializing the full CDS position list.

    Parameters
    ----------
    enst_id : str
        Ensembl transcript ID (without version suffix).
    residue_positions : list of int
        1-based protein residue positions to map.
    cds_lookup : dict
        Mapping of ENST -> sorted list of (chrom, start, end, strand, phase).
    cache : dict, optional
        Shared TranscriptCDS cache for cross-call reuse.

    Returns
    -------
    list of (int, CodonMapping or None)
        Tuples of (residue_pos, mapping_or_None).
    """
    tcds = _get_transcript_cds(enst_id, cds_lookup, cache)
    if tcds is None:
        return [(p, None) for p in residue_positions]
    return [(p, tcds.map_residue(p)) for p in residue_positions]


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
