"""Single-pass Ensembl GTF -> per-gene structural summary.

Deliberately separate from ``hvantk.algorithms.ptm.mapper.parse_ensembl_gtf``: that one
keys MANE by gene SYMBOL and builds per-transcript exon lists for residue-to-genomic
coordinate mapping. This one keys by ``gene_id`` and emits one row per gene. Pure Python
and pandas -- no Hail -- so it stays in the fast test suite.

The representative coding transcript is MANE Select when the gene has one, otherwise the
transcript with the longest summed CDS. Length is preferred over "first seen" because GTF
transcript order is not stable across releases.
"""
from __future__ import annotations

import gzip
import logging
import re
from collections import defaultdict

import pandas as pd

logger = logging.getLogger(__name__)

_RE_GENE_ID = re.compile(r'gene_id "([^"]+)"')
_RE_TX_ID = re.compile(r'transcript_id "([^"]+)"')
_RE_BIOTYPE = re.compile(r'gene_biotype "([^"]+)"')

COLUMNS = [
    "gene_id",
    "gene_biotype",
    "mane_select",
    "cds_transcript",
    "cds_length",
    "n_coding_exons",
    "n_transcripts",
]


def parse_gtf_structure(gtf_path: str) -> pd.DataFrame:
    """Summarise an Ensembl GTF to one row per gene.

    Parameters
    ----------
    gtf_path : str
        Path to an Ensembl GTF, gzipped or plain.

    Returns
    -------
    pandas.DataFrame
        Columns as in ``COLUMNS``, sorted by ``gene_id`` with a reset index.
    """
    tx_of_gene: dict[str, set[str]] = defaultdict(set)
    cds_bp: dict[str, int] = defaultdict(int)
    cds_n: dict[str, int] = defaultdict(int)
    mane_of_gene: dict[str, str] = {}
    biotype_of_gene: dict[str, str] = {}

    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if len(fields) < 9:
                continue
            feature, attrs = fields[2], fields[8]
            if feature != "transcript" and feature != "CDS":
                continue

            m_gene = _RE_GENE_ID.search(attrs)
            m_tx = _RE_TX_ID.search(attrs)
            if not m_gene or not m_tx:
                continue
            gene_id = m_gene.group(1).split(".")[0]
            tx_id = m_tx.group(1).split(".")[0]

            if feature == "transcript":
                tx_of_gene[gene_id].add(tx_id)
                if gene_id not in biotype_of_gene:
                    m_bt = _RE_BIOTYPE.search(attrs)
                    if m_bt:
                        biotype_of_gene[gene_id] = m_bt.group(1)
                if 'tag "MANE_Select"' in attrs:
                    mane_of_gene[gene_id] = tx_id
            else:  # CDS -- GTF coordinates are 1-based inclusive on both ends
                cds_bp[tx_id] += int(fields[4]) - int(fields[3]) + 1
                cds_n[tx_id] += 1

    rows = []
    for gene_id, transcripts in tx_of_gene.items():
        coding = [t for t in transcripts if cds_bp[t] > 0]
        representative = mane_of_gene.get(gene_id)
        if representative is None or representative not in coding:
            # Tie-break on transcript ID, not just length: `coding` derives from a set, so a
            # bare max() on length alone resolves ties by set iteration order, which follows
            # Python's per-process string-hash seed. That makes the output differ between runs
            # on identical input.
            representative = (
                max(coding, key=lambda t: (cds_bp[t], t)) if coding else None
            )
        rows.append(
            {
                "gene_id": gene_id,
                "gene_biotype": biotype_of_gene.get(gene_id, ""),
                "mane_select": mane_of_gene.get(gene_id, ""),
                "cds_transcript": representative or "",
                "cds_length": cds_bp[representative] if representative else 0,
                "n_coding_exons": cds_n[representative] if representative else 0,
                "n_transcripts": len(transcripts),
            }
        )

    df = pd.DataFrame(rows, columns=COLUMNS)
    logger.info("Parsed %d genes from %s", len(df), gtf_path)
    return df.sort_values("gene_id").reset_index(drop=True)
