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
import io
import logging
import os
import re
import shutil
import urllib.parse
import zipfile
from collections import defaultdict
from dataclasses import dataclass
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

# Phospho modification masses in PeptideAtlas bracket notation.
# These are total modified residue masses (amino acid + HPO3):
#   Phosphoserine:   S(87.03) + HPO3(79.97) ≈ 167.0
#   Phosphothreonine: T(101.05) + HPO3(79.97) ≈ 181.0
#   Phosphotyrosine:  Y(163.06) + HPO3(79.97) ≈ 243.0
_PHOSPHO_BRACKET_MASS_BY_AA = {
    "S": 167.0,
    "T": 181.0,
    "Y": 243.0,
}
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


def _iter_tsv_from_zip(zf: zipfile.ZipFile, table_name: str):
    """Iterate over rows of a TSV table inside a zip without loading all into memory."""
    member = _find_table_in_zip(zf, table_name)
    if member is None:
        logger.warning("Table %s not found in zip archive", table_name)
        return

    with zf.open(member) as f:
        text = io.TextIOWrapper(f, encoding="utf-8")
        reader = csv.DictReader(text, delimiter="\t")
        yield from reader


def _extract_phospho_offsets(modified_sequence: str) -> List[tuple]:
    """Extract 0-based offsets and amino acids of phosphorylated S/T/Y.

    PeptideAtlas uses two notations for modifications:
      - Text: ``S[Phospho]``, ``T[Phospho]``, ``Y[Phospho]``
      - Numeric mass: ``S[167]``, ``T[181]``, ``Y[243]``

    N-terminal labels like ``[TMT6plex]-`` are skipped (no preceding residue).

    Returns a list of (0-based_offset, amino_acid) tuples for phospho sites.
    """
    offsets = []
    residue_index = -1  # will be incremented to 0 for first AA
    last_aa = None  # track the preceding amino acid
    i = 0
    seq = modified_sequence
    while i < len(seq):
        ch = seq[i]
        if ch == "[":
            try:
                end = seq.index("]", i)
            except ValueError:
                break  # Malformed — unclosed bracket
            bracket_content = seq[i + 1 : end]
            is_phospho = (
                bracket_content == "Phospho"
                or bracket_content.startswith("Phospho:")
            )
            # Also handle numeric mass notation for S/T/Y.
            if not is_phospho:
                try:
                    mass = float(bracket_content)
                    expected_mass = _PHOSPHO_BRACKET_MASS_BY_AA.get(last_aa)
                    if (
                        expected_mass is not None
                        and abs(mass - expected_mass) <= _PHOSPHO_BRACKET_TOLERANCE
                    ):
                        is_phospho = True
                except ValueError:
                    pass

            # Only accept phospho on S, T, Y residues
            if is_phospho and residue_index >= 0 and last_aa in _PHOSPHO_AA_DESC:
                offsets.append((residue_index, last_aa))

            i = end + 1
        elif ch.isalpha() and ch.isupper():
            residue_index += 1
            last_aa = ch
            i += 1
        else:
            # Skip lowercase, digits, dashes, etc.
            i += 1
    return offsets


def parse_peptideatlas_zip(zip_path: str) -> List[dict]:  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    """Parse a PeptideAtlas TSV zip and extract phospho sites.

    Joins biosequence, peptide_instance, peptide_mapping, and
    modified_peptide_instance tables. Filters out DECOY_ and CONTAM_
    prefixed accessions. Aggregates observation counts for the same
    protein site across multiple peptide instances.

    Uses a streaming approach: builds lightweight lookup dicts for
    smaller tables (biosequence, peptide_instance), then streams the
    large tables (peptide_mapping, modified_peptide_instance) without
    loading them fully into memory.

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
        required_tables = [
            "biosequence.tsv",
            "peptide_instance.tsv",
            "peptide_mapping.tsv",
            "modified_peptide_instance.tsv",
        ]
        missing = [t for t in required_tables if _find_table_in_zip(zf, t) is None]
        if missing:
            raise FileNotFoundError(
                f"Missing required table(s) in {zip_path}: {', '.join(missing)}. "
                f"Available files: {zf.namelist()}"
            )

        # Step 0: Index canonical proteins from protein_identification table
        # presence_level_id=1 is "canonical" in PeptideAtlas
        logger.info("Indexing canonical proteins...")
        canonical_bs_ids: set = set()
        for row in _iter_tsv_from_zip(zf, "protein_identification.tsv"):
            if row.get("presence_level_id") == "1":
                canonical_bs_ids.add(row["biosequence_id"])
        logger.info("Found %d canonical proteins", len(canonical_bs_ids))

        # Step 1: Index biosequences — only keep canonical, non-DECOY proteins
        # Store minimal info: {biosequence_id: (accession, gene_name, seq_len)}
        logger.info("Indexing biosequences...")
        bioseq_by_id: Dict[str, tuple] = {}
        n_bio = 0
        for row in _iter_tsv_from_zip(zf, "biosequence.tsv"):
            n_bio += 1
            bs_id = row["biosequence_id"]
            acc = row.get("biosequence_accession", "")
            if acc.startswith("DECOY_") or acc.startswith("CONTAM_"):
                continue
            if canonical_bs_ids and bs_id not in canonical_bs_ids:
                continue
            bioseq_by_id[bs_id] = (
                acc,
                row.get("biosequence_gene_name", ""),
                str(len(row.get("biosequence_seq", ""))),
            )
        logger.info("Indexed %d biosequences (%d canonical)", n_bio, len(bioseq_by_id))

        # Step 2: Index peptide_instance (381K rows)
        # Store: {peptide_instance_id: n_observations}
        logger.info("Indexing peptide instances...")
        pi_obs: Dict[str, int] = {}
        for row in _iter_tsv_from_zip(zf, "peptide_instance.tsv"):
            pi_obs[row["peptide_instance_id"]] = int(row.get("n_observations", "0"))
        logger.info("Indexed %d peptide instances", len(pi_obs))

        # Step 3: Stream peptide_mapping (15M rows) — build a compact index
        # Store: {peptide_instance_id: [(biosequence_id, start_in_biosequence), ...]}
        logger.info("Streaming peptide mappings...")
        mapping_by_pi: Dict[str, List[tuple]] = defaultdict(list)
        n_mappings = 0
        for row in _iter_tsv_from_zip(zf, "peptide_mapping.tsv"):
            n_mappings += 1
            pi_id = row["peptide_instance_id"]
            bs_id = row.get("matched_biosequence_id", "")
            if bs_id in bioseq_by_id:  # Only keep mappings to non-decoy proteins
                start = int(row["start_in_biosequence"])
                mapping_by_pi[pi_id].append((bs_id, start))
            if n_mappings % 5_000_000 == 0:
                logger.info("  ...processed %dM peptide mappings", n_mappings // 1_000_000)
        logger.info("Processed %d peptide mappings, %d with valid proteins",
                     n_mappings, len(mapping_by_pi))

        # Step 4: Stream modified_peptide_instance (3M rows) — extract phospho sites
        logger.info("Extracting phospho sites from modified peptide instances...")
        site_obs: Dict[tuple, int] = defaultdict(int)
        site_info: Dict[tuple, dict] = {}
        n_mpi = 0
        n_with_phospho = 0

        for mp in _iter_tsv_from_zip(zf, "modified_peptide_instance.tsv"):
            n_mpi += 1
            pi_id = mp["peptide_instance_id"]
            if pi_id not in pi_obs:
                continue

            mod_seq = mp.get("modified_peptide_sequence", "")
            phospho_sites = _extract_phospho_offsets(mod_seq)
            if not phospho_sites:
                continue

            n_with_phospho += 1
            n_obs = pi_obs[pi_id]
            mappings = mapping_by_pi.get(pi_id, [])

            for bs_id, start in mappings:
                bs = bioseq_by_id[bs_id]
                acc, gene_name, seq_len = bs

                for offset, aa_char in phospho_sites:
                    site_pos = start + offset
                    key = (acc, site_pos)
                    site_obs[key] += n_obs

                    if key not in site_info:
                        site_info[key] = {
                            "accession": acc,
                            "gene_symbol": gene_name,
                            "position": site_pos,
                            "amino_acid": aa_char,
                            "description": _PHOSPHO_AA_DESC.get(aa_char, "Phosphorylation"),
                            "ensembl_xrefs": "",
                            "sequence_length": seq_len,
                        }

            if n_mpi % 1_000_000 == 0:
                logger.info("  ...processed %dM modified peptides (%d with phospho)",
                             n_mpi // 1_000_000, n_with_phospho)

    # Build final list
    sites = []
    for key, info in site_info.items():
        site = dict(info)
        site["n_observations"] = site_obs[key]
        sites.append(site)

    logger.info("Parsed %d distinct phospho sites from %d modified peptides "
                 "(%d with phospho)", len(sites), n_mpi, n_with_phospho)
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

    with open(output_path, "w", newline="", encoding="utf-8") as f:
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

    @staticmethod
    def _validate_zip_url(url: str) -> None:
        """Validate that download URL uses HTTPS and points to peptideatlas.org."""
        parsed = urllib.parse.urlparse(url)
        if parsed.scheme != "https":
            raise ValueError(
                "Only HTTPS download URLs are allowed (got %s): %s"
                % (parsed.scheme or "<empty>", url)
            )
        host = parsed.hostname or ""
        if host != "peptideatlas.org" and not host.endswith(".peptideatlas.org"):
            raise ValueError(
                "Unexpected download host for PeptideAtlas URL "
                "(expected peptideatlas.org or *.peptideatlas.org): %s"
                % url
            )

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

        Notes
        -----
        The generated URL follows:
        ``{PEPTIDEATLAS_PHOSPHO_BASE_URL}/{build_date}/atlas_build_{build_id}.tsv.zip``.
        """
        zip_url = (
            f"{PEPTIDEATLAS_PHOSPHO_BASE_URL}/{build_date}/"
            f"atlas_build_{build_id}.tsv.zip"
        )
        return cls(
            build_date=build_date,
            build_id=build_id,
            zip_url=zip_url,
        )

    def download(self, output_dir: str, overwrite: bool = False) -> str:
        """Download PeptideAtlas phospho build and produce intermediate TSV.

        Downloads the TSV zip, parses phospho sites with observation counts
        from canonical proteins, and writes an intermediate TSV compatible
        with the PTM pipeline mapper.

        Parameters
        ----------
        output_dir : str
            Directory to save files (zip + intermediate TSV).
        overwrite : bool
            If True, re-download and re-parse even if files exist.

        Returns
        -------
        str
            Path to the intermediate TSV file.
        """
        os.makedirs(output_dir, exist_ok=True)
        zip_filename = f"atlas_build_{self.build_id}.tsv.zip"
        zip_path = os.path.join(output_dir, zip_filename)
        tsv_filename = f"peptideatlas-phospho-{self.build_date}-{self.build_id}.tsv"
        tsv_path = os.path.join(output_dir, tsv_filename)

        if os.path.exists(tsv_path) and not overwrite:
            logger.info("Intermediate TSV already exists: %s", tsv_path)
            return tsv_path

        # Download zip if needed
        if not os.path.exists(zip_path) or overwrite:
            self._validate_zip_url(self.zip_url)
            logger.info("Downloading %s -> %s", self.zip_url, zip_path)
            import requests

            with requests.get(self.zip_url, stream=True, timeout=600) as resp:
                resp.raise_for_status()
                with open(zip_path, "wb") as fout:
                    shutil.copyfileobj(resp.raw, fout)
            logger.info("Download complete: %s", zip_path)
        else:
            logger.info("Using cached zip: %s", zip_path)

        # Parse and write intermediate TSV
        logger.info("Parsing phospho sites from zip...")
        sites = parse_peptideatlas_zip(zip_path)
        write_intermediate_tsv(sites, tsv_path)
        return tsv_path

    def get_metadata(self) -> dict:
        """Return metadata about this dataset."""
        return {
            "source": "PeptideAtlas Human Phospho Build",
            "build_date": self.build_date,
            "build_id": self.build_id,
            "zip_url": self.zip_url,
        }
