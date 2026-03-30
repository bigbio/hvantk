# PeptideAtlas Phospho Handler Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add PeptideAtlas human phospho build as a data source in the PTM pipeline, downloading the TSV dump, parsing phospho sites with observation counts, and mapping them to genomic coordinates in the 13-column stream format.

**Architecture:** New dataset class (`PeptideAtlasPhosphoDataset`) downloads and parses the PeptideAtlas TSV zip, aggregates phospho sites per protein/position with counts, and outputs an intermediate TSV. The existing GTF-based mapper translates protein positions to genomic codons. The output column schema extends from 12 to 13 columns (adding `n_observations`). Both UniProt and PeptideAtlas sources feed into the same Hail Table builder.

**Tech Stack:** Python 3.10+, Click (CLI), csv module, zipfile, pytest, Hail (table builder)

**Spec:** `docs/superpowers/specs/2026-03-30-peptideatlas-phospho-handler-design.md`

---

## File Structure

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `hvantk/datasets/peptideatlas_phospho_datasets.py` | Dataset class: download zip, parse tables, aggregate phospho sites, write intermediate TSV |
| Create | `hvantk/commands/peptideatlas_phospho_downloader.py` | CLI command `hvantk download peptideatlas-phospho` |
| Create | `hvantk/tests/test_peptideatlas_phospho.py` | Unit tests for dataset class parsing and aggregation |
| Modify | `hvantk/ptm/constants.py` | Add PeptideAtlas constants, extend `PTM_OUTPUT_COLUMNS` to 13 cols |
| Modify | `hvantk/ptm/pipeline.py` | Update `map_ptm_sites()` for 13 cols, add PeptideAtlas path to `PTMBuildConfig` and `ptm_build_pipeline()` |
| Modify | `hvantk/commands/ptm_cli.py` | Add `--peptideatlas-tsv` option to `ptm build` |
| Modify | `hvantk/commands/download_cli.py` | Register new downloader |
| Modify | `hvantk/tables/table_builders.py` | Add `n_observations` to import types and transform |
| Modify | `hvantk/tests/test_ptm.py` | Update round-trip test for 13-column output |

---

### Task 1: Extend PTM_OUTPUT_COLUMNS and Add PeptideAtlas Constants

**Files:**
- Modify: `hvantk/ptm/constants.py:1-59`

- [ ] **Step 1: Add PeptideAtlas constants and extend PTM_OUTPUT_COLUMNS**

In `hvantk/ptm/constants.py`, add PeptideAtlas URL constants after the Ensembl GTF block (after line 15), and update `PTM_OUTPUT_COLUMNS` to include `n_observations`:

```python
# PeptideAtlas Phospho Build
PEPTIDEATLAS_PHOSPHO_BASE_URL = "http://www.peptideatlas.org/builds"
PEPTIDEATLAS_LATEST_BUILD_DATE = "202512"
PEPTIDEATLAS_LATEST_BUILD_ID = "606"
```

Update `PTM_OUTPUT_COLUMNS` (currently lines 42-55) to:

```python
PTM_OUTPUT_COLUMNS = [
    "chrom",
    "codon_start",
    "codon_end",
    "strand",
    "uniprot_id",
    "gene_symbol",
    "residue_pos",
    "amino_acid",
    "ptm_type",
    "ptm_category",
    "source_db",
    "evidence_type",
    "n_observations",
]
```

- [ ] **Step 2: Run existing tests to confirm no breakage**

Run: `pytest hvantk/tests/test_ptm.py -v`

Expected: `test_map_ptm_sites_roundtrip` will FAIL because the writer now expects 13 columns but the `writerow()` in `pipeline.py` only produces 12. This is expected — we fix it in Task 2.

- [ ] **Step 3: Commit**

```bash
git add hvantk/ptm/constants.py
git commit -m "feat: add PeptideAtlas constants and n_observations to PTM_OUTPUT_COLUMNS"
```

---

### Task 2: Update Pipeline to Write 13-Column Output

**Files:**
- Modify: `hvantk/ptm/pipeline.py:150-265`
- Modify: `hvantk/tests/test_ptm.py:130-199`

- [ ] **Step 1: Update `map_ptm_sites()` writer to include `n_observations`**

In `hvantk/ptm/pipeline.py`, update the `writer.writerow()` call inside `flush_protein()` (around line 225) to pass through `source_db`, `evidence_type`, and `n_observations` from the input record, with defaults for UniProt:

```python
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
```

- [ ] **Step 2: Update the round-trip test for 13 columns**

In `hvantk/tests/test_ptm.py`, update `test_map_ptm_sites_roundtrip` to verify the new columns. After the existing assertion `assert ser315["source_db"] == "UniProt"` (line 196), add:

```python
    assert ser315["evidence_type"] == "curated"
    assert ser315["n_observations"] == "0"
```

- [ ] **Step 3: Run tests to verify they pass**

Run: `pytest hvantk/tests/test_ptm.py -v`
Expected: ALL PASS (including `test_map_ptm_sites_roundtrip`)

- [ ] **Step 4: Commit**

```bash
git add hvantk/ptm/pipeline.py hvantk/tests/test_ptm.py
git commit -m "feat: update PTM pipeline to write 13-column output with n_observations"
```

---

### Task 3: Update Hail Table Builder for n_observations

**Files:**
- Modify: `hvantk/tables/table_builders.py:1357-1464`

- [ ] **Step 1: Add `n_observations` to import types and transform**

In `hvantk/tables/table_builders.py`, update `create_ptm_sites_tb`:

1. Update the docstring (line 1371) to mention the 13th column:
```python
    with columns: chrom, codon_start, codon_end, strand, uniprot_id,
    gene_symbol, residue_pos, amino_acid, ptm_type, ptm_category,
    source_db, evidence_type, n_observations.
```

2. Add `n_observations` to the import types dict (line 1454-1458):
```python
            types={
                "codon_start": hl.tstr,
                "codon_end": hl.tstr,
                "residue_pos": hl.tstr,
                "n_observations": hl.tstr,
            },
```

3. Add `n_observations` cast in the transform function, after the existing `hl.int32` casts (line 1422-1427):
```python
        # Cast numeric fields
        ht = ht.annotate(
            codon_start=hl.int32(ht.codon_start),
            codon_end=hl.int32(ht.codon_end),
            residue_pos=hl.int32(ht.residue_pos),
            n_observations=hl.int32(ht.n_observations),
        )
```

- [ ] **Step 2: Commit**

```bash
git add hvantk/tables/table_builders.py
git commit -m "feat: add n_observations field to PTM sites Hail Table builder"
```

---

### Task 4: Create PeptideAtlas Dataset Class — Test First

**Files:**
- Create: `hvantk/tests/test_peptideatlas_phospho.py`
- Create: `hvantk/datasets/peptideatlas_phospho_datasets.py`

- [ ] **Step 1: Write failing tests for PeptideAtlas parsing**

Create `hvantk/tests/test_peptideatlas_phospho.py`:

```python
"""Tests for hvantk.datasets.peptideatlas_phospho_datasets module.

Covers:
1. Phospho site extraction from mock PeptideAtlas TSV tables
2. Observation count aggregation (multiple peptides → same site)
3. Dataset class construction and metadata
4. Intermediate TSV output format
"""

import csv
import os
import zipfile

import pytest


# ---------- Fixtures ----------


def _write_tsv(path, fieldnames, rows):
    """Helper to write a TSV file with given columns and rows."""
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


@pytest.fixture
def mock_pa_zip(tmp_path):
    """Create a mock PeptideAtlas TSV zip with minimal phospho data.

    Simulates: TP53 protein with two phospho sites (Ser315, Ser6),
    where Ser315 is observed by two distinct peptides (counts should sum).
    """
    tables_dir = tmp_path / "tables"
    tables_dir.mkdir()

    # biosequence table — maps biosequence_id to UniProt accession
    _write_tsv(
        tables_dir / "biosequence.tsv",
        ["biosequence_id", "biosequence_name", "biosequence_accession",
         "biosequence_gene_name", "biosequence_seq", "organism_id"],
        [
            {
                "biosequence_id": "100",
                "biosequence_name": "TP53_HUMAN",
                "biosequence_accession": "P04637",
                "biosequence_gene_name": "TP53",
                "biosequence_seq": "M" * 393,
                "organism_id": "9606",
            },
            {
                "biosequence_id": "200",
                "biosequence_name": "DECOY_FOO",
                "biosequence_accession": "DECOY_Q99999",
                "biosequence_gene_name": "FOO",
                "biosequence_seq": "M" * 100,
                "organism_id": "9606",
            },
        ],
    )

    # peptide_instance table — peptide instances with observation counts
    _write_tsv(
        tables_dir / "peptide_instance.tsv",
        ["peptide_instance_id", "peptide_id", "n_observations", "n_samples"],
        [
            {"peptide_instance_id": "1", "peptide_id": "10", "n_observations": "50", "n_samples": "5"},
            {"peptide_instance_id": "2", "peptide_id": "20", "n_observations": "30", "n_samples": "3"},
            {"peptide_instance_id": "3", "peptide_id": "30", "n_observations": "10", "n_samples": "1"},
        ],
    )

    # peptide_mapping table — maps peptide_instance to protein position
    _write_tsv(
        tables_dir / "peptide_mapping.tsv",
        ["peptide_instance_id", "matched_biosequence_id", "start_in_biosequence", "end_in_biosequence"],
        [
            {"peptide_instance_id": "1", "matched_biosequence_id": "100", "start_in_biosequence": "310", "end_in_biosequence": "320"},
            {"peptide_instance_id": "2", "matched_biosequence_id": "100", "start_in_biosequence": "312", "end_in_biosequence": "325"},
            {"peptide_instance_id": "3", "matched_biosequence_id": "100", "start_in_biosequence": "1", "end_in_biosequence": "15"},
        ],
    )

    # modified_peptide_instance table — phospho modifications within peptides
    # modification_offset is 0-based within the peptide sequence
    _write_tsv(
        tables_dir / "modified_peptide_instance.tsv",
        ["modified_peptide_instance_id", "peptide_instance_id",
         "modified_peptide_sequence", "modification_mass"],
        [
            # Peptide 1: phos at offset 5 within peptide starting at pos 310 → site 315
            {
                "modified_peptide_instance_id": "1001",
                "peptide_instance_id": "1",
                "modified_peptide_sequence": "AAAAAS[167]AAAAA",
                "modification_mass": "79.9663",
            },
            # Peptide 2: phos at offset 3 within peptide starting at pos 312 → site 315 (same site!)
            {
                "modified_peptide_instance_id": "1002",
                "peptide_instance_id": "2",
                "modified_peptide_sequence": "AAS[167]AAAAAAAAA",
                "modification_mass": "79.9663",
            },
            # Peptide 3: phos at offset 5 within peptide starting at pos 1 → site 6
            {
                "modified_peptide_instance_id": "1003",
                "peptide_instance_id": "3",
                "modified_peptide_sequence": "AAAAAS[167]AAAAAAAAA",
                "modification_mass": "79.9663",
            },
        ],
    )

    # Create zip
    zip_path = tmp_path / "atlas_build_test.tsv.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        for tsv_file in tables_dir.glob("*.tsv"):
            zf.write(tsv_file, tsv_file.name)

    return zip_path


# ---------- Test 1: Phospho site extraction ----------


def test_parse_phospho_sites(mock_pa_zip, tmp_path):
    """Phospho sites are correctly extracted with positions and counts."""
    from hvantk.datasets.peptideatlas_phospho_datasets import parse_peptideatlas_zip

    sites = parse_peptideatlas_zip(str(mock_pa_zip))

    # Should find 2 distinct sites on TP53: position 315 and position 6
    assert len(sites) == 2

    by_pos = {s["position"]: s for s in sites}
    assert 315 in by_pos
    assert 6 in by_pos

    # Ser315 was observed by peptide 1 (50 obs) and peptide 2 (30 obs) = 80 total
    assert by_pos[315]["n_observations"] == 80
    assert by_pos[315]["accession"] == "P04637"
    assert by_pos[315]["gene_symbol"] == "TP53"

    # Ser6 was observed by peptide 3 only (10 obs)
    assert by_pos[6]["n_observations"] == 10


# ---------- Test 2: DECOY filtering ----------


def test_decoy_sequences_filtered(mock_pa_zip):
    """DECOY and contaminant sequences are excluded."""
    from hvantk.datasets.peptideatlas_phospho_datasets import parse_peptideatlas_zip

    sites = parse_peptideatlas_zip(str(mock_pa_zip))
    accessions = {s["accession"] for s in sites}
    assert all(not acc.startswith("DECOY_") for acc in accessions)


# ---------- Test 3: Dataset class ----------


def test_dataset_from_latest():
    """from_latest() creates a dataset with current build defaults."""
    from hvantk.datasets.peptideatlas_phospho_datasets import PeptideAtlasPhosphoDataset

    dataset = PeptideAtlasPhosphoDataset.from_latest()
    assert dataset.build_date == "202512"
    assert dataset.build_id == "606"
    meta = dataset.get_metadata()
    assert meta["source"] == "PeptideAtlas Human Phospho Build"


def test_dataset_from_build():
    """from_build() creates a dataset with specified build parameters."""
    from hvantk.datasets.peptideatlas_phospho_datasets import PeptideAtlasPhosphoDataset

    dataset = PeptideAtlasPhosphoDataset.from_build("202204", "500")
    assert dataset.build_date == "202204"
    assert dataset.build_id == "500"
    assert "202204" in dataset.zip_url


# ---------- Test 4: Intermediate TSV output ----------


def test_write_intermediate_tsv(mock_pa_zip, tmp_path):
    """write_intermediate_tsv produces correctly formatted TSV."""
    from hvantk.datasets.peptideatlas_phospho_datasets import (
        parse_peptideatlas_zip,
        write_intermediate_tsv,
    )

    sites = parse_peptideatlas_zip(str(mock_pa_zip))
    output_path = str(tmp_path / "peptideatlas_phospho.tsv")
    write_intermediate_tsv(sites, output_path)

    with open(output_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 2
    expected_cols = {
        "accession", "gene_symbol", "position", "description",
        "amino_acid", "ensembl_xrefs", "sequence_length",
        "n_observations", "source_db", "evidence_type",
    }
    assert set(reader.fieldnames) == expected_cols

    row315 = next(r for r in rows if r["position"] == "315")
    assert row315["source_db"] == "PeptideAtlas"
    assert row315["evidence_type"] == "mass_spectrometry"
    assert row315["n_observations"] == "80"
    assert row315["description"] == "Phosphoserine"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest hvantk/tests/test_peptideatlas_phospho.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'hvantk.datasets.peptideatlas_phospho_datasets'`

- [ ] **Step 3: Implement the dataset class**

Create `hvantk/datasets/peptideatlas_phospho_datasets.py`:

```python
"""
PeptideAtlas Human Phospho Build dataset handling.

Downloads and parses the PeptideAtlas phospho build TSV dump to extract
phosphorylation site positions with observation counts. Produces an
intermediate TSV compatible with the hvantk PTM pipeline mapper.

Example usage:
    dataset = PeptideAtlasPhosphoDataset.from_latest()
    dataset.download("/data/peptideatlas/")

    # Or parse a pre-downloaded zip:
    sites = parse_peptideatlas_zip("/data/peptideatlas/atlas_build_606.tsv.zip")
    write_intermediate_tsv(sites, "/data/peptideatlas/phospho_sites.tsv")
"""

import csv
import logging
import os
import re
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

# Intermediate TSV columns (compatible with UniProt PTM TSV + extras)
_INTERMEDIATE_COLUMNS = [
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

# Phospho modification mass (monoisotopic, +/- tolerance for matching)
_PHOSPHO_MASS = 79.9663
_PHOSPHO_MASS_TOL = 0.01

# Regex to find modification brackets in modified peptide sequences
# Matches patterns like S[167], T[181], Y[243]
_MOD_BRACKET_RE = re.compile(r"([A-Z])\[([0-9.]+)\]")

# Phospho-modified amino acids and their PTM descriptions
_PHOSPHO_AA_DESC = {
    "S": "Phosphoserine",
    "T": "Phosphothreonine",
    "Y": "Phosphotyrosine",
}


def _find_table_in_zip(zf: zipfile.ZipFile, table_name: str) -> Optional[str]:
    """Find a TSV file in the zip matching the given table name.

    Searches for exact match first, then substring match. Returns
    the zip member name or None.
    """
    names = zf.namelist()
    # Exact match (e.g., "biosequence.tsv")
    for name in names:
        basename = os.path.basename(name).lower()
        if basename == f"{table_name}.tsv":
            return name
    # Substring match (e.g., "atlas_build_606_biosequence.tsv")
    for name in names:
        basename = os.path.basename(name).lower()
        if table_name in basename and basename.endswith(".tsv"):
            return name
    return None


def _read_table(zf: zipfile.ZipFile, member_name: str) -> List[dict]:
    """Read a TSV table from inside a zip, returning list of row dicts."""
    with zf.open(member_name) as f:
        import io
        text = io.TextIOWrapper(f, encoding="utf-8")
        reader = csv.DictReader(text, delimiter="\t")
        return list(reader)


def _extract_phospho_offsets(modified_sequence: str) -> List[int]:
    """Extract 0-based offsets of phosphorylation sites within a peptide.

    Parses the PeptideAtlas modified peptide sequence notation where
    modifications are encoded as amino_acid[mass]. Phosphorylation is
    identified by mass ~79.966 Da on S, T, or Y.

    Args:
        modified_sequence: e.g., "AAAAAS[167]AAAAA"

    Returns:
        List of 0-based positions within the *unmodified* peptide
        where phospho modifications occur.
    """
    offsets = []
    # Track position in the unmodified sequence
    clean_pos = 0
    i = 0
    seq = modified_sequence

    while i < len(seq):
        ch = seq[i]
        if ch.isalpha() and ch.isupper():
            # Check if next char is '[' (modification)
            if i + 1 < len(seq) and seq[i + 1] == "[":
                # Find closing bracket
                bracket_end = seq.index("]", i + 2)
                mass_str = seq[i + 2:bracket_end]
                try:
                    mass = float(mass_str)
                except ValueError:
                    mass = 0.0

                # Check if this is a phospho modification on S/T/Y
                if ch in _PHOSPHO_AA_DESC and abs(mass - 167.0) < 1.0:
                    offsets.append(clean_pos)

                clean_pos += 1
                i = bracket_end + 1
            else:
                clean_pos += 1
                i += 1
        else:
            i += 1

    return offsets


def parse_peptideatlas_zip(zip_path: str) -> List[dict]:
    """Parse a PeptideAtlas TSV dump zip and extract phospho sites.

    Reads the relevant database tables from the zip, joins them to
    determine phospho site positions on proteins, and aggregates
    observation counts per (accession, position).

    Args:
        zip_path: Path to the PeptideAtlas TSV zip file.

    Returns:
        List of dicts with keys: accession, gene_symbol, position,
        amino_acid, sequence_length, n_observations.

    Raises:
        FileNotFoundError: If required tables are missing from the zip.
    """
    with zipfile.ZipFile(zip_path, "r") as zf:
        # Discover tables
        available = zf.namelist()
        logger.info(f"PeptideAtlas zip contains {len(available)} files: {available}")

        bio_name = _find_table_in_zip(zf, "biosequence")
        pi_name = _find_table_in_zip(zf, "peptide_instance")
        pm_name = _find_table_in_zip(zf, "peptide_mapping")
        mpi_name = _find_table_in_zip(zf, "modified_peptide_instance")

        missing = []
        if not bio_name:
            missing.append("biosequence")
        if not pi_name:
            missing.append("peptide_instance")
        if not pm_name:
            missing.append("peptide_mapping")
        if not mpi_name:
            missing.append("modified_peptide_instance")

        if missing:
            raise FileNotFoundError(
                f"Required tables not found in zip: {missing}. "
                f"Available files: {available}"
            )

        logger.info(f"Reading tables: {bio_name}, {pi_name}, {pm_name}, {mpi_name}")

        # 1. Parse biosequence → {biosequence_id: {accession, gene, seq_len}}
        bio_rows = _read_table(zf, bio_name)
        biosequences = {}
        for row in bio_rows:
            acc = row.get("biosequence_accession", "")
            # Skip DECOY and contaminant sequences
            if acc.startswith("DECOY_") or acc.startswith("CONTAM_"):
                continue
            bid = row.get("biosequence_id", "")
            biosequences[bid] = {
                "accession": acc,
                "gene_symbol": row.get("biosequence_gene_name", ""),
                "sequence_length": str(len(row.get("biosequence_seq", ""))),
            }
        logger.info(f"Parsed {len(biosequences)} biosequences (excluding decoys)")

        # 2. Parse peptide_instance → {peptide_instance_id: n_observations}
        pi_rows = _read_table(zf, pi_name)
        instance_obs = {}
        for row in pi_rows:
            piid = row.get("peptide_instance_id", "")
            n_obs = int(row.get("n_observations", "0"))
            instance_obs[piid] = n_obs
        logger.info(f"Parsed {len(instance_obs)} peptide instances")

        # 3. Parse peptide_mapping → {peptide_instance_id: [(bio_id, start)]}
        pm_rows = _read_table(zf, pm_name)
        instance_mappings = defaultdict(list)
        for row in pm_rows:
            piid = row.get("peptide_instance_id", "")
            bio_id = row.get("matched_biosequence_id", "")
            start = int(row.get("start_in_biosequence", "0"))
            instance_mappings[piid].append((bio_id, start))
        logger.info(f"Parsed {len(pm_rows)} peptide mappings")

        # 4. Parse modified_peptide_instance → phospho offsets per peptide_instance
        mpi_rows = _read_table(zf, mpi_name)
        instance_phospho_offsets = defaultdict(list)
        for row in mpi_rows:
            piid = row.get("peptide_instance_id", "")
            mod_seq = row.get("modified_peptide_sequence", "")
            offsets = _extract_phospho_offsets(mod_seq)
            if offsets:
                instance_phospho_offsets[piid].extend(offsets)
        logger.info(
            f"Parsed {len(mpi_rows)} modified peptide instances, "
            f"{len(instance_phospho_offsets)} with phospho sites"
        )

    # 5. Join: for each modified peptide instance with phospho sites,
    #    map to protein position and aggregate counts
    # Key: (accession, position) → total n_observations
    site_counts = defaultdict(int)
    site_info = {}

    for piid, offsets in instance_phospho_offsets.items():
        n_obs = instance_obs.get(piid, 0)
        mappings = instance_mappings.get(piid, [])

        for bio_id, pep_start in mappings:
            bio = biosequences.get(bio_id)
            if bio is None:
                continue  # Decoy or unknown

            acc = bio["accession"]
            for offset in offsets:
                site_pos = pep_start + offset
                key = (acc, site_pos)
                site_counts[key] += n_obs
                if key not in site_info:
                    site_info[key] = {
                        "accession": acc,
                        "gene_symbol": bio["gene_symbol"],
                        "sequence_length": bio["sequence_length"],
                    }

    # 6. Build output list
    sites = []
    for (acc, pos), count in sorted(site_counts.items()):
        info = site_info[(acc, pos)]
        sites.append(
            {
                "accession": info["accession"],
                "gene_symbol": info["gene_symbol"],
                "position": pos,
                "amino_acid": "",  # Not available from PeptideAtlas tables directly
                "sequence_length": info["sequence_length"],
                "n_observations": count,
            }
        )

    logger.info(f"Aggregated {len(sites)} distinct phospho sites")
    return sites


def write_intermediate_tsv(sites: List[dict], output_path: str) -> str:
    """Write parsed phospho sites to an intermediate TSV.

    The output is compatible with the PTM pipeline mapper input format,
    with additional columns for source_db, evidence_type, and n_observations.

    Args:
        sites: List of site dicts from parse_peptideatlas_zip().
        output_path: Path to write the TSV file.

    Returns:
        Path to the written file.
    """
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=_INTERMEDIATE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()

        for site in sites:
            aa = site.get("amino_acid", "")
            desc = _PHOSPHO_AA_DESC.get(aa, "Phosphoserine")
            writer.writerow(
                {
                    "accession": site["accession"],
                    "gene_symbol": site["gene_symbol"],
                    "position": str(site["position"]),
                    "description": desc,
                    "amino_acid": aa,
                    "ensembl_xrefs": "",  # Resolved via gene_mane fallback
                    "sequence_length": site["sequence_length"],
                    "n_observations": str(site["n_observations"]),
                    "source_db": "PeptideAtlas",
                    "evidence_type": "mass_spectrometry",
                }
            )

    logger.info(f"Wrote {len(sites)} phospho sites to {output_path}")
    return output_path


@dataclass
class PeptideAtlasPhosphoDataset:
    """Represents a PeptideAtlas human phospho build dataset.

    Attributes:
        build_date: Build date in YYYYMM format (e.g., "202512").
        build_id: Build ID number (e.g., "606").
        zip_url: URL to the TSV zip file.
        file_name: Name of the output intermediate TSV file.
    """

    build_date: str
    build_id: str
    zip_url: str
    file_name: str

    @classmethod
    def from_build(cls, build_date: str, build_id: str) -> "PeptideAtlasPhosphoDataset":
        """Create a dataset reference for a specific build.

        Args:
            build_date: Build date in YYYYMM format.
            build_id: Build ID number.

        Returns:
            PeptideAtlasPhosphoDataset instance.
        """
        zip_url = (
            f"{PEPTIDEATLAS_PHOSPHO_BASE_URL}/{build_date}"
            f"/atlas_build_{build_id}.tsv.zip"
        )
        file_name = f"peptideatlas-phospho-{build_date}.tsv"

        return cls(
            build_date=build_date,
            build_id=build_id,
            zip_url=zip_url,
            file_name=file_name,
        )

    @classmethod
    def from_latest(cls) -> "PeptideAtlasPhosphoDataset":
        """Create a dataset reference for the latest available build.

        Returns:
            PeptideAtlasPhosphoDataset instance for the latest build.
        """
        return cls.from_build(
            PEPTIDEATLAS_LATEST_BUILD_DATE,
            PEPTIDEATLAS_LATEST_BUILD_ID,
        )

    def download(self, output_dir: str, overwrite: bool = False) -> str:
        """Download and parse the PeptideAtlas phospho build.

        Downloads the TSV zip, parses phospho sites with counts,
        and writes an intermediate TSV for the PTM pipeline.

        Args:
            output_dir: Directory to save files.
            overwrite: If True, overwrite existing files.

        Returns:
            Path to the intermediate TSV file.

        Raises:
            FileExistsError: If file exists and overwrite is False.
            RuntimeError: If download or parsing fails.
        """
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, self.file_name)

        if os.path.exists(output_path) and not overwrite:
            raise FileExistsError(
                f"File already exists: {output_path}. Use overwrite=True to replace."
            )

        zip_path = os.path.join(output_dir, os.path.basename(self.zip_url))

        try:
            # Download zip
            logger.info(f"Downloading PeptideAtlas phospho build from {self.zip_url}")
            from urllib.request import urlretrieve
            urlretrieve(self.zip_url, zip_path)
            logger.info(f"Downloaded: {zip_path}")

            # Parse and write intermediate TSV
            sites = parse_peptideatlas_zip(zip_path)
            write_intermediate_tsv(sites, output_path)

            return output_path

        except Exception as e:
            if os.path.exists(output_path):
                os.remove(output_path)
            raise RuntimeError(
                f"Failed to download/parse PeptideAtlas phospho data: {e}"
            ) from e

    def get_metadata(self) -> Dict[str, str]:
        """Get metadata about this dataset.

        Returns:
            Dictionary with dataset metadata.
        """
        return {
            "source": "PeptideAtlas Human Phospho Build",
            "build_date": self.build_date,
            "build_id": self.build_id,
            "zip_url": self.zip_url,
            "file_name": self.file_name,
            "description": (
                "Phosphorylation sites from the PeptideAtlas human phospho "
                "build, with observation counts (PSMs) per site. Sites are "
                "aggregated across all peptides covering each position."
            ),
        }

    def __str__(self) -> str:
        return f"PeptideAtlasPhosphoDataset(build={self.build_date}, id={self.build_id})"
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest hvantk/tests/test_peptideatlas_phospho.py -v`
Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add hvantk/datasets/peptideatlas_phospho_datasets.py hvantk/tests/test_peptideatlas_phospho.py
git commit -m "feat: add PeptideAtlas phospho dataset class with parsing and tests"
```

---

### Task 5: Create CLI Downloader and Register It

**Files:**
- Create: `hvantk/commands/peptideatlas_phospho_downloader.py`
- Modify: `hvantk/commands/download_cli.py:1-24`

- [ ] **Step 1: Create the downloader CLI command**

Create `hvantk/commands/peptideatlas_phospho_downloader.py`:

```python
"""
CLI command to download PeptideAtlas human phospho build data.

Example:
    hvantk download peptideatlas-phospho --output-dir /data/peptideatlas/
"""

import logging

import click

logger = logging.getLogger(__name__)


@click.command(name="peptideatlas-phospho")
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    required=True,
    help="Directory to save the PeptideAtlas phospho TSV file",
)
@click.option(
    "--build-date",
    type=str,
    default=None,
    help="Build date in YYYYMM format (default: latest build)",
)
@click.option(
    "--build-id",
    type=str,
    default=None,
    help="Build ID number (default: latest build)",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing file",
)
@click.pass_context
def peptideatlas_phospho_downloader(ctx, output_dir, build_date, build_id, overwrite):
    """Download human phospho site data from PeptideAtlas.

    Downloads the PeptideAtlas phospho build TSV dump, parses phospho
    site positions with observation counts, and writes an intermediate
    TSV compatible with the PTM pipeline.
    """
    try:
        from hvantk.datasets.peptideatlas_phospho_datasets import PeptideAtlasPhosphoDataset

        if build_date and build_id:
            dataset = PeptideAtlasPhosphoDataset.from_build(build_date, build_id)
        else:
            dataset = PeptideAtlasPhosphoDataset.from_latest()

        click.echo(f"Downloading PeptideAtlas phospho data ({dataset})...")
        output_path = dataset.download(output_dir, overwrite=overwrite)
        click.echo(f"Downloaded to: {output_path}")

    except Exception as e:
        logger.exception(f"Download failed: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
```

- [ ] **Step 2: Register the downloader in download_cli.py**

In `hvantk/commands/download_cli.py`, add the import (after line 10) and registration (after line 24):

Import:
```python
from hvantk.commands.peptideatlas_phospho_downloader import peptideatlas_phospho_downloader
```

Registration:
```python
download_group.add_command(peptideatlas_phospho_downloader, "peptideatlas-phospho")
```

- [ ] **Step 3: Test CLI registration**

Run: `python -m hvantk.hvantk download --help`
Expected: Output includes `peptideatlas-phospho` in the list of subcommands.

- [ ] **Step 4: Commit**

```bash
git add hvantk/commands/peptideatlas_phospho_downloader.py hvantk/commands/download_cli.py
git commit -m "feat: add peptideatlas-phospho download CLI command"
```

---

### Task 6: Integrate PeptideAtlas into PTM Build Pipeline

**Files:**
- Modify: `hvantk/ptm/pipeline.py:40-336`
- Modify: `hvantk/commands/ptm_cli.py:54-145`

- [ ] **Step 1: Add `peptideatlas_tsv` to PTMBuildConfig**

In `hvantk/ptm/pipeline.py`, add to the `PTMBuildConfig` dataclass (after the `ptm_tsv` field at line 53):

```python
    peptideatlas_tsv: Optional[str] = None
```

Update `validate()` to check the new field (after line 72):

```python
        if self.peptideatlas_tsv and not os.path.exists(self.peptideatlas_tsv):
            errors.append(f"PeptideAtlas TSV file not found: {self.peptideatlas_tsv}")
```

- [ ] **Step 2: Update `ptm_build_pipeline()` to handle PeptideAtlas source**

In `hvantk/ptm/pipeline.py`, after the existing `map_ptm_sites` call (line 314), before `# Step 5: Build Hail Table` (line 317), add PeptideAtlas mapping and concatenation:

```python
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
```

Also add `PTM_OUTPUT_COLUMNS` to the imports at the top of the file (line 28) if not already imported — it's already imported at line 28.

- [ ] **Step 3: Add `--peptideatlas-tsv` option to CLI build command**

In `hvantk/commands/ptm_cli.py`, add a new option to the `ptm_build` command (after the `--ptm-tsv` option block, around line 80):

```python
@click.option(
    "--peptideatlas-tsv",
    type=click.Path(exists=True),
    default=None,
    help="Path to PeptideAtlas phospho intermediate TSV (from 'hvantk download peptideatlas-phospho')",
)
```

Update the function signature (line 93) to include the new parameter:

```python
def ptm_build(
    ctx, output_dir, output_ht, gtf_path, ptm_tsv, peptideatlas_tsv, flanking_codons, overwrite
):
```

Pass it through to the config (around line 115):

```python
        config = PTMBuildConfig(
            output_dir=output_dir,
            output_ht=output_ht,
            gtf_path=gtf_path,
            ptm_tsv=ptm_tsv,
            peptideatlas_tsv=peptideatlas_tsv,
            flanking_codons=flanking_codons,
            overwrite=overwrite,
        )
```

- [ ] **Step 4: Run all tests**

Run: `pytest hvantk/tests/test_ptm.py hvantk/tests/test_peptideatlas_phospho.py -v`
Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add hvantk/ptm/pipeline.py hvantk/commands/ptm_cli.py
git commit -m "feat: integrate PeptideAtlas phospho source into PTM build pipeline"
```

---

### Task 7: Final Verification

- [ ] **Step 1: Run full test suite**

Run: `pytest hvantk/tests/ -v --tb=short`
Expected: ALL PASS

- [ ] **Step 2: Verify CLI help renders correctly**

Run: `python -m hvantk.hvantk ptm build --help`
Expected: Output shows `--peptideatlas-tsv` option alongside existing options.

Run: `python -m hvantk.hvantk download peptideatlas-phospho --help`
Expected: Output shows `--output-dir`, `--build-date`, `--build-id`, `--overwrite` options.

- [ ] **Step 3: Commit any final fixes**

If any tests failed or help output was wrong, fix and commit.
