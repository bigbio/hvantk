"""
PeptideAtlas Phospho dataset handling.

This module provides functions and a dataset class for downloading and parsing
human phospho-proteome data from PeptideAtlas builds. PeptideAtlas aggregates
mass-spectrometry observations of phosphorylated peptides across hundreds of
experiments.

Example usage:
    dataset = PeptideAtlasPhosphoDataset.from_latest()
    dataset.download("/data/peptideatlas")

    sites = parse_peptideatlas_zip("/data/peptideatlas/atlas_build_606.tsv.zip")
    write_intermediate_tsv(sites, "/data/peptideatlas/peptideatlas_phospho.tsv")
"""

import csv
import logging
import os
import re
import zipfile
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from hvantk.ptm.constants import (
    PEPTIDEATLAS_PHOSPHO_BASE_URL,
    PEPTIDEATLAS_LATEST_BUILD_DATE,
    PEPTIDEATLAS_LATEST_BUILD_ID,
)

logger = logging.getLogger(__name__)

# Phospho amino acid descriptions
_PHOSPHO_AA_DESC = {
    "S": "Phosphoserine",
    "T": "Phosphothreonine",
    "Y": "Phosphotyrosine",
}

# Phospho modification mass in PeptideAtlas bracket notation (e.g. S[167])
_PHOSPHO_BRACKET_MASS = 167.0
_PHOSPHO_BRACKET_TOLERANCE = 1.0

# Output TSV column order for intermediate file
_TSV_COLUMNS = [
    "accession",
    "gene_symbol",
    "position",
    "description",
    "amino_acid",
    "ensembl_xrefs",
    "sequence_length",
    "n_observations",
    "source_db",
    "evidence_type",
]

# Pattern to match bracket modifications like [167] or [167.0]
_BRACKET_MOD_PATTERN = re.compile(r"\[([0-9.]+)\]")


def _find_table_in_zip(zf: zipfile.ZipFile, table_name: str) -> Optional[str]:
    """Find a table file in a PeptideAtlas zip by name.

    Tries exact match first (e.g. ``biosequence.tsv``), then substring match
    (e.g. ``atlas_build_606_biosequence.tsv``).
    """
    names = zf.namelist()
    # Exact match
    if table_name in names:
        return table_name
    # Substring match
    for name in names:
        if table_name in name:
            return name
    return None


def _read_tsv_from_zip(zf: zipfile.ZipFile, table_name: str) -> List[Dict[str, str]]:
    """Read a TSV table from within a zip file, returning list of row dicts."""
    member = _find_table_in_zip(zf, table_name)
    if member is None:
        logger.warning("Table %s not found in zip archive", table_name)
        return []
    with zf.open(member) as f:
        text = f.read().decode("utf-8")
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    return list(reader)


def _extract_phospho_offsets(modified_sequence: str) -> List[int]:
    """Extract 0-based offsets of phosphorylated S/T/Y within a peptide.

    PeptideAtlas notation: ``AAAAAS[167]AAAAA`` where ``[167]`` indicates
    a phospho modification on the preceding residue.

    Returns a list of 0-based residue positions within the peptide.
    """
    offsets = []
    residue_index = -1  # will be incremented to 0 for first AA
    i = 0
    seq = modified_sequence
    while i < len(seq):
        ch = seq[i]
        if ch == "[":
            # Extract the mass inside brackets
            end = seq.index("]", i)
            mass_str = seq[i + 1 : end]
            try:
                mass = float(mass_str)
            except ValueError:
                i = end + 1
                continue
            # Check if this is a phospho modification
            if abs(mass - _PHOSPHO_BRACKET_MASS) <= _PHOSPHO_BRACKET_TOLERANCE:
                # The modification applies to the preceding residue
                if residue_index >= 0:
                    preceding_aa = None
                    # Walk back to find the preceding AA character
                    # residue_index is the 0-based index of the last AA seen
                    # We need the actual character
                    aa_count = 0
                    for c in modified_sequence:
                        if c == "[":
                            break_at = modified_sequence.index("]", modified_sequence.index("[")) + 1
                        if c.isalpha():
                            if aa_count == residue_index:
                                preceding_aa = c
                                break
                            aa_count += 1
                    if preceding_aa and preceding_aa.upper() in _PHOSPHO_AA_DESC:
                        offsets.append(residue_index)
            i = end + 1
        elif ch.isalpha():
            residue_index += 1
            i += 1
        else:
            i += 1
    return offsets


def parse_peptideatlas_zip(zip_path: str) -> List[dict]:
    """Parse a PeptideAtlas TSV zip and extract phospho sites.

    Joins biosequence, peptide_instance, peptide_mapping, and
    modified_peptide_instance tables. Filters out DECOY_ and CONTAM_
    prefixed accessions. Aggregates observation counts for the same
    protein site across multiple peptide instances.

    Parameters
    ----------
    zip_path : str
        Path to the PeptideAtlas ``*.tsv.zip`` archive.

    Returns
    -------
    list of dict
        Each dict has keys: accession, gene_symbol, position, amino_acid,
        description, ensembl_xrefs, sequence_length, n_observations.
    """
    with zipfile.ZipFile(zip_path, "r") as zf:
        biosequences = _read_tsv_from_zip(zf, "biosequence.tsv")
        peptide_instances = _read_tsv_from_zip(zf, "peptide_instance.tsv")
        peptide_mappings = _read_tsv_from_zip(zf, "peptide_mapping.tsv")
        modified_peptides = _read_tsv_from_zip(zf, "modified_peptide_instance.tsv")

    # Index biosequences by id, filtering DECOY_ and CONTAM_
    bioseq_by_id: Dict[str, dict] = {}
    for bs in biosequences:
        acc = bs.get("biosequence_accession", "")
        if acc.startswith("DECOY_") or acc.startswith("CONTAM_"):
            continue
        bioseq_by_id[bs["biosequence_id"]] = bs

    # Index peptide instances by id
    pi_by_id: Dict[str, dict] = {}
    for pi in peptide_instances:
        pi_by_id[pi["peptide_instance_id"]] = pi

    # Index peptide mappings by peptide_instance_id
    mapping_by_pi: Dict[str, List[dict]] = defaultdict(list)
    for pm in peptide_mappings:
        mapping_by_pi[pm["peptide_instance_id"]].append(pm)

    # Aggregate phospho sites: key = (accession, position) -> total observations
    site_obs: Dict[tuple, int] = defaultdict(int)
    site_info: Dict[tuple, dict] = {}

    for mp in modified_peptides:
        pi_id = mp["peptide_instance_id"]
        if pi_id not in pi_by_id:
            continue

        offsets = _extract_phospho_offsets(mp["modified_peptide_sequence"])
        if not offsets:
            continue

        n_obs = int(pi_by_id[pi_id].get("n_observations", 0))
        mappings = mapping_by_pi.get(pi_id, [])

        for mapping in mappings:
            bs_id = mapping["matched_biosequence_id"]
            if bs_id not in bioseq_by_id:
                continue

            bs = bioseq_by_id[bs_id]
            start = int(mapping["start_in_biosequence"])
            acc = bs["biosequence_accession"]

            for offset in offsets:
                site_pos = start + offset
                key = (acc, site_pos)
                site_obs[key] += n_obs

                if key not in site_info:
                    # Determine the amino acid at this offset
                    mod_seq = mp["modified_peptide_sequence"]
                    aa_count = 0
                    aa_char = "S"  # default
                    for c in mod_seq:
                        if c.isalpha():
                            if aa_count == offset:
                                aa_char = c.upper()
                                break
                            aa_count += 1

                    site_info[key] = {
                        "accession": acc,
                        "gene_symbol": bs.get("biosequence_gene_name", ""),
                        "position": site_pos,
                        "amino_acid": aa_char,
                        "description": _PHOSPHO_AA_DESC.get(aa_char, "Phosphorylation"),
                        "ensembl_xrefs": "",
                        "sequence_length": len(bs.get("biosequence_seq", "")),
                    }

    # Build final list
    sites = []
    for key, info in site_info.items():
        site = dict(info)
        site["n_observations"] = site_obs[key]
        sites.append(site)

    logger.info("Parsed %d phospho sites from %s", len(sites), zip_path)
    return sites


def write_intermediate_tsv(sites: List[dict], output_path: str) -> str:
    """Write phospho sites to an intermediate TSV file.

    Parameters
    ----------
    sites : list of dict
        Output from :func:`parse_peptideatlas_zip`.
    output_path : str
        Destination file path.

    Returns
    -------
    str
        The output file path.
    """
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=_TSV_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for site in sites:
            row = dict(site)
            row["source_db"] = "PeptideAtlas"
            row["evidence_type"] = "mass_spectrometry"
            writer.writerow(row)

    logger.info("Wrote %d phospho sites to %s", len(sites), output_path)
    return output_path


@dataclass
class PeptideAtlasPhosphoDataset:
    """Represents a PeptideAtlas Human Phospho build.

    Use :meth:`from_latest` or :meth:`from_build` to construct instances.
    """

    build_date: str
    build_id: str
    zip_url: str

    @classmethod
    def from_latest(cls) -> "PeptideAtlasPhosphoDataset":
        """Create a dataset pointing to the latest known build."""
        return cls.from_build(
            PEPTIDEATLAS_LATEST_BUILD_DATE,
            PEPTIDEATLAS_LATEST_BUILD_ID,
        )

    @classmethod
    def from_build(cls, build_date: str, build_id: str) -> "PeptideAtlasPhosphoDataset":
        """Create a dataset for a specific PeptideAtlas build.

        Parameters
        ----------
        build_date : str
            Build date string, e.g. ``"202512"``.
        build_id : str
            Build numeric ID, e.g. ``"606"``.
        """
        zip_url = (
            f"{PEPTIDEATLAS_PHOSPHO_BASE_URL}/phospho/{build_date}/"
            f"atlas_build_{build_id}.tsv.zip"
        )
        return cls(
            build_date=build_date,
            build_id=build_id,
            zip_url=zip_url,
        )

    def download(self, output_dir: str, overwrite: bool = False) -> str:
        """Download the PeptideAtlas zip to output_dir.

        Parameters
        ----------
        output_dir : str
            Directory to save the zip file.
        overwrite : bool
            If True, re-download even if the file exists.

        Returns
        -------
        str
            Path to the downloaded zip file.
        """
        import urllib.request

        os.makedirs(output_dir, exist_ok=True)
        filename = f"atlas_build_{self.build_id}.tsv.zip"
        output_path = os.path.join(output_dir, filename)

        if os.path.exists(output_path) and not overwrite:
            logger.info("File already exists: %s", output_path)
            return output_path

        logger.info("Downloading %s -> %s", self.zip_url, output_path)
        urllib.request.urlretrieve(self.zip_url, output_path)
        logger.info("Download complete: %s", output_path)
        return output_path

    def get_metadata(self) -> dict:
        """Return metadata about this dataset."""
        return {
            "source": "PeptideAtlas Human Phospho Build",
            "build_date": self.build_date,
            "build_id": self.build_id,
            "zip_url": self.zip_url,
        }
