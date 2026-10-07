"""Tests for hvantk.skills.peptideatlas.phospho.shared.datasets module.

Covers:
1. Phospho site extraction from mock PeptideAtlas TSV tables
2. Observation count aggregation (multiple peptides -> same site)
3. Dataset class construction and metadata
4. Intermediate TSV output format
5. _extract_phospho_offsets parser (text/numeric notation, edge cases)
6. Canonical protein filtering via protein_identification table
"""

import csv
import os
import zipfile
from unittest import mock

import pytest


# ---------- Fixtures ----------


def _write_tsv(path, fieldnames, rows):
    """Helper to write a TSV file with given columns and rows."""
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _rewrite_zip_table(zip_path, tmp_path, table_name, replacements):
    rewritten = tmp_path / f"rewritten-{zip_path.name}"
    with zipfile.ZipFile(zip_path) as source, zipfile.ZipFile(rewritten, "w") as target:
        for name in source.namelist():
            content = source.read(name)
            if name == table_name:
                text = content.decode("utf-8")
                for old, new in replacements:
                    text = text.replace(old, new)
                content = text.encode("utf-8")
            target.writestr(name, content)
    return rewritten


@pytest.fixture
def mock_pa_zip(tmp_path):
    """Create a mock PeptideAtlas TSV zip with minimal phospho data.

    Simulates: TP53 protein with two phospho sites (Ser315, Ser6),
    where Ser315 is observed by two distinct peptides (counts should sum).

    Also maps peptide_instance 4 (phospho, modified_peptide_instance 1004) onto
    the DECOY_Q99999 biosequence (id 200), so the DECOY_ accession-prefix filter
    actually has a phospho site to remove. Without this mapping, DECOY_Q99999
    carries no peptide at all and `test_decoy_sequences_filtered`'s assertion
    passes vacuously regardless of whether the prefix filter exists.
    """
    tables_dir = tmp_path / "tables"
    tables_dir.mkdir()

    # biosequence table
    _write_tsv(
        tables_dir / "biosequence.tsv",
        [
            "biosequence_id",
            "biosequence_name",
            "biosequence_accession",
            "biosequence_gene_name",
            "biosequence_seq",
            "organism_id",
        ],
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

    # peptide_instance table
    _write_tsv(
        tables_dir / "peptide_instance.tsv",
        ["peptide_instance_id", "peptide_id", "n_observations", "n_samples"],
        [
            {
                "peptide_instance_id": "1",
                "peptide_id": "10",
                "n_observations": "50",
                "n_samples": "5",
            },
            {
                "peptide_instance_id": "2",
                "peptide_id": "20",
                "n_observations": "30",
                "n_samples": "3",
            },
            {
                "peptide_instance_id": "3",
                "peptide_id": "30",
                "n_observations": "10",
                "n_samples": "1",
            },
            {
                "peptide_instance_id": "4",
                "peptide_id": "40",
                "n_observations": "15",
                "n_samples": "2",
            },
        ],
    )

    # peptide_mapping table
    _write_tsv(
        tables_dir / "peptide_mapping.tsv",
        [
            "peptide_instance_id",
            "matched_biosequence_id",
            "start_in_biosequence",
            "end_in_biosequence",
        ],
        [
            {
                "peptide_instance_id": "1",
                "matched_biosequence_id": "100",
                "start_in_biosequence": "310",
                "end_in_biosequence": "320",
            },
            {
                "peptide_instance_id": "2",
                "matched_biosequence_id": "100",
                "start_in_biosequence": "313",
                "end_in_biosequence": "325",
            },
            {
                "peptide_instance_id": "3",
                "matched_biosequence_id": "100",
                "start_in_biosequence": "1",
                "end_in_biosequence": "15",
            },
            {
                "peptide_instance_id": "4",
                "matched_biosequence_id": "200",
                "start_in_biosequence": "5",
                "end_in_biosequence": "15",
            },
        ],
    )

    # modified_peptide_instance table
    _write_tsv(
        tables_dir / "modified_peptide_instance.tsv",
        [
            "modified_peptide_instance_id",
            "peptide_instance_id",
            "modified_peptide_sequence",
            "modification_mass",
            "n_observations",
        ],
        [
            {
                "modified_peptide_instance_id": "1001",
                "peptide_instance_id": "1",
                "modified_peptide_sequence": "AAAAAS[167]AAAAA",
                "modification_mass": "79.9663",
                "n_observations": "50",
            },
            {
                "modified_peptide_instance_id": "1002",
                "peptide_instance_id": "2",
                "modified_peptide_sequence": "AAS[167]AAAAAAAAA",
                "modification_mass": "79.9663",
                "n_observations": "30",
            },
            {
                "modified_peptide_instance_id": "1003",
                "peptide_instance_id": "3",
                "modified_peptide_sequence": "AAAAAS[167]AAAAAAAAA",
                "modification_mass": "79.9663",
                "n_observations": "10",
            },
            {
                "modified_peptide_instance_id": "1004",
                "peptide_instance_id": "4",
                "modified_peptide_sequence": "AAAAS[Phospho]AAAAA",
                "modification_mass": "79.9663",
                "n_observations": "15",
            },
        ],
    )

    zip_path = tmp_path / "atlas_build_test.tsv.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        for tsv_file in tables_dir.glob("*.tsv"):
            zf.write(tsv_file, tsv_file.name)

    return zip_path


# ---------- Test 1: Phospho site extraction ----------


def test_parse_phospho_sites(mock_pa_zip, tmp_path):
    """Phospho sites are correctly extracted with positions and counts."""
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        parse_peptideatlas_zip,
    )

    sites = parse_peptideatlas_zip(str(mock_pa_zip))

    assert len(sites) == 2

    by_pos = {s["position"]: s for s in sites}
    assert 315 in by_pos
    assert 6 in by_pos

    # Ser315: peptide 1 (50) + peptide 2 (30) = 80
    assert by_pos[315]["n_observations"] == 80
    assert by_pos[315]["accession"] == "P04637"
    assert by_pos[315]["gene_symbol"] == "TP53"

    assert by_pos[6]["n_observations"] == 10


def _two_form_zip(tmp_path, form_counts, unmodified_count="80", unmapped_count=None):
    """Peptide 1, seen 100 times, mapped to P04637: a phospho form per count in
    ``form_counts`` plus its unmodified form (``unmodified_count``). With
    ``unmapped_count``, peptide 2 (no mapping to a kept protein) adds one phospho
    form with that count."""
    tables = tmp_path / "two_forms"
    tables.mkdir()
    _write_tsv(
        tables / "biosequence.tsv",
        [
            "biosequence_id",
            "biosequence_accession",
            "biosequence_gene_name",
            "biosequence_seq",
        ],
        [
            {
                "biosequence_id": "100",
                "biosequence_accession": "P04637",
                "biosequence_gene_name": "TP53",
                "biosequence_seq": "M" * 393,
            }
        ],
    )
    _write_tsv(
        tables / "peptide_instance.tsv",
        ["peptide_instance_id", "n_observations"],
        [
            {"peptide_instance_id": "1", "n_observations": "100"},
            {"peptide_instance_id": "2", "n_observations": "5"},
        ],
    )
    _write_tsv(
        tables / "peptide_mapping.tsv",
        ["peptide_instance_id", "matched_biosequence_id", "start_in_biosequence"],
        [
            {
                "peptide_instance_id": "1",
                "matched_biosequence_id": "100",
                "start_in_biosequence": "310",
            }
        ],
    )
    _write_tsv(
        tables / "modified_peptide_instance.tsv",
        [
            "modified_peptide_instance_id",
            "peptide_instance_id",
            "modified_peptide_sequence",
            "peptide_charge",
            "n_observations",
        ],
        [
            {
                "modified_peptide_instance_id": str(1000 + charge),
                "peptide_instance_id": "1",
                "modified_peptide_sequence": "AAAAAS[Phospho]AAAAA",
                "peptide_charge": str(charge),
                "n_observations": count,
            }
            for charge, count in enumerate(form_counts, start=2)
        ]
        + [
            {
                "modified_peptide_instance_id": "1100",
                "peptide_instance_id": "1",
                "modified_peptide_sequence": "AAAAASAAAAA",
                "peptide_charge": "2",
                "n_observations": unmodified_count,
            }
        ]
        + (
            []
            if unmapped_count is None
            else [
                {
                    "modified_peptide_instance_id": "2001",
                    "peptide_instance_id": "2",
                    "modified_peptide_sequence": "AAS[Phospho]AA",
                    "peptide_charge": "2",
                    "n_observations": unmapped_count,
                }
            ]
        ),
    )
    zip_path = tmp_path / "two_forms.tsv.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        for tsv_file in tables.glob("*.tsv"):
            zf.write(tsv_file, tsv_file.name)
    return str(zip_path)


def test_site_counts_sum_each_forms_own_observations(tmp_path):
    """#425: phospho forms seen 12 and 8 times give the site 20: not the peptide's
    100 once per form, and not its 100 over all forms (the unmodified form, 80,
    does not carry the site)."""
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        parse_peptideatlas_zip,
    )

    sites = parse_peptideatlas_zip(_two_form_zip(tmp_path, ["12", "8"]))
    assert [(s["position"], s["n_observations"]) for s in sites] == [(315, 20)]


def test_form_without_a_count_fails(tmp_path):
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        parse_peptideatlas_zip,
    )

    with pytest.raises(ValueError, match="has no integer n_observations"):
        parse_peptideatlas_zip(_two_form_zip(tmp_path, ["12", "\\N"]))


@pytest.mark.parametrize(
    "kwargs", [{"unmodified_count": "\\N"}, {"unmapped_count": "\\N"}]
)
def test_forms_that_add_nothing_are_not_checked(tmp_path, kwargs):
    """A count is read only from a phospho form that maps to a kept protein: an
    unmodified form, or a form of a peptide with no kept mapping, may lack one."""
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        parse_peptideatlas_zip,
    )

    sites = parse_peptideatlas_zip(_two_form_zip(tmp_path, ["12", "8"], **kwargs))
    assert [(s["position"], s["n_observations"]) for s in sites] == [(315, 20)]


# ---------- Test 2: DECOY filtering ----------


def test_decoy_sequences_filtered(mock_pa_zip):
    """DECOY and contaminant sequences are excluded.

    `mock_pa_zip` maps a real phospho peptide onto the DECOY_Q99999 biosequence,
    so this exercises the DECOY_ accession-prefix filter itself (dropping a site
    that would otherwise be produced), not just the absence of DECOY input.
    """
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        parse_peptideatlas_zip,
    )

    sites = parse_peptideatlas_zip(str(mock_pa_zip))
    accessions = {s["accession"] for s in sites}
    assert all(not acc.startswith("DECOY_") for acc in accessions)
    assert len(sites) == 2, "the DECOY-mapped site must not appear in the output"


# ---------- Test 3: Dataset class ----------


def test_dataset_from_latest():
    """from_latest() creates a dataset with current build defaults."""
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        PeptideAtlasPhosphoDataset,
    )

    dataset = PeptideAtlasPhosphoDataset.from_latest()
    assert dataset.build_date == "202512"
    assert dataset.build_id == "606"
    meta = dataset.get_metadata()
    assert meta["source"] == "PeptideAtlas Human Phospho Build"


def test_dataset_from_build():
    """from_build() creates a dataset with specified build parameters."""
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        PeptideAtlasPhosphoDataset,
    )

    dataset = PeptideAtlasPhosphoDataset.from_build("202204", "500")
    assert dataset.build_date == "202204"
    assert dataset.build_id == "500"
    assert "202204" in dataset.zip_url
    assert "/phospho/202204/" in dataset.zip_url


# ---------- Test 4: Intermediate TSV output ----------


def test_write_intermediate_tsv(mock_pa_zip, tmp_path):
    """write_intermediate_tsv produces correctly formatted TSV."""
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
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
    }
    assert set(reader.fieldnames) == expected_cols

    row315 = next(r for r in rows if r["position"] == "315")
    assert row315["source_db"] == "PeptideAtlas"
    assert row315["evidence_type"] == "mass_spectrometry"
    assert row315["n_observations"] == "80"
    assert row315["description"] == "Phosphoserine"


# ---------- Test 5: _extract_phospho_offsets ----------


@pytest.mark.parametrize(
    "mod_seq, expected",
    [
        # Text notation
        ("DSY[Phospho]VGDEAQSK", [(2, "Y")]),
        ("DTC[Carbamidomethyl]YS[Phospho]PK", [(4, "S")]),
        ("TVT[Phospho]PLNQVANPNSAIFGGARPR", [(2, "T")]),
        # Numeric mass notation (modified residue masses)
        ("AAAAAS[167]AAAAA", [(5, "S")]),
        ("AAAAT[181]AAAAA", [(4, "T")]),
        ("AAAAY[243]AAAAA", [(4, "Y")]),
        ("AS[+79.966]K", [(1, "S")]),
        ("AS[+79.97]K", [(1, "S")]),
        ("AS[+80.5]K", []),
        ("AT[UNIMOD:21]K", [(1, "T")]),
        # Sulfation (79.9568 Da) is close enough to phospho (79.96633 Da) that
        # a loose tolerance wrongly accepted it; must stay rejected.
        ("AY[+79.957]K", []),
        # N-terminal label skipped (no preceding residue)
        ("[TMT6plex]-DSY[Phospho]VGDEAQSK", [(2, "Y")]),
        ("[iTRAQ4plex]-AAST[Phospho]R", [(3, "T")]),
        # Non-S/T/Y residue — should NOT be phospho
        ("A[167]R", []),
        ("DVA[Phospho]TPLNQR", []),
        # Non-phospho modification on S/T/Y — should NOT be phospho (build 606
        # has none, so the fixture cannot carry one; these unit cases can)
        ("AS[Oxidation]K", []),
        ("AAS[Acetyl]AAT[Phospho]K", [(5, "T")]),
        ("AAS[129]AA", []),
        # Multiple modifications on same peptide
        ("AS[Phospho]AT[Phospho]Y", [(1, "S"), (3, "T")]),
        ("S[Phospho]T[Phospho]Y[Phospho]", [(0, "S"), (1, "T"), (2, "Y")]),
        # Malformed bracket — should not crash
        ("AS[", []),
        # No modifications
        ("PEPTIDER", []),
        # Empty string
        ("", []),
    ],
)
def test_extract_phospho_offsets(mod_seq, expected):
    """_extract_phospho_offsets handles all PeptideAtlas notations."""
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        _extract_phospho_offsets,
    )

    result = _extract_phospho_offsets(mod_seq)
    assert result == expected


# ---------- Test 6: Missing tables ----------


def test_parse_raises_when_required_tables_missing(tmp_path):
    """Missing required tables should raise FileNotFoundError."""
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        parse_peptideatlas_zip,
    )

    zip_path = tmp_path / "missing_tables.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr(
            "biosequence.tsv", "biosequence_id\tbiosequence_accession\n1\tP12345\n"
        )

    with pytest.raises(FileNotFoundError):
        parse_peptideatlas_zip(str(zip_path))


# ---------- Test 7: Dataset download ----------


def test_dataset_download_returns_intermediate_tsv(tmp_path, mock_pa_zip):
    """download() returns generated intermediate TSV path."""
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        PeptideAtlasPhosphoDataset,
    )

    dataset = PeptideAtlasPhosphoDataset.from_build("202512", "606")
    output_dir = tmp_path / "out"

    class _FakeResp:
        def __init__(self, src_path):
            self._src = open(src_path, "rb")
            self.raw = self._src

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self._src.close()

        def raise_for_status(self):
            return None

    with mock.patch(
        "requests.get",
        side_effect=lambda url, **kwargs: _FakeResp(mock_pa_zip),
    ):
        tsv_path = dataset.download(str(output_dir), overwrite=True)

    assert tsv_path.endswith("peptideatlas-phospho-202512-606.tsv")
    assert os.path.exists(tsv_path)


# ---------- Test 8: Canonical protein filtering ----------


def test_canonical_protein_filtering(tmp_path):
    """Only canonical proteins (presence_level_id=1) are included."""
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        parse_peptideatlas_zip,
    )

    tables_dir = tmp_path / "tables"
    tables_dir.mkdir()

    # biosequence: two proteins
    _write_tsv(
        tables_dir / "biosequence.tsv",
        [
            "biosequence_id",
            "biosequence_name",
            "biosequence_accession",
            "biosequence_gene_name",
            "biosequence_seq",
            "organism_id",
        ],
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
                "biosequence_id": "300",
                "biosequence_name": "FOO_HUMAN",
                "biosequence_accession": "Q99999",
                "biosequence_gene_name": "FOO",
                "biosequence_seq": "M" * 200,
                "organism_id": "9606",
            },
        ],
    )

    # protein_identification: only biosequence 100 is canonical
    _write_tsv(
        tables_dir / "protein_identification.tsv",
        ["biosequence_id", "presence_level_id", "protein_identification_id"],
        [
            {
                "biosequence_id": "100",
                "presence_level_id": "1",
                "protein_identification_id": "1",
            },
            {
                "biosequence_id": "300",
                "presence_level_id": "3",
                "protein_identification_id": "2",
            },
        ],
    )

    # peptide_instance
    _write_tsv(
        tables_dir / "peptide_instance.tsv",
        ["peptide_instance_id", "peptide_id", "n_observations", "n_samples"],
        [
            {
                "peptide_instance_id": "1",
                "peptide_id": "10",
                "n_observations": "50",
                "n_samples": "5",
            },
        ],
    )

    # peptide_mapping: peptide maps to BOTH proteins
    _write_tsv(
        tables_dir / "peptide_mapping.tsv",
        [
            "peptide_instance_id",
            "matched_biosequence_id",
            "start_in_biosequence",
            "end_in_biosequence",
        ],
        [
            {
                "peptide_instance_id": "1",
                "matched_biosequence_id": "100",
                "start_in_biosequence": "310",
                "end_in_biosequence": "320",
            },
            {
                "peptide_instance_id": "1",
                "matched_biosequence_id": "300",
                "start_in_biosequence": "50",
                "end_in_biosequence": "60",
            },
        ],
    )

    # modified_peptide_instance: phospho at offset 5
    _write_tsv(
        tables_dir / "modified_peptide_instance.tsv",
        [
            "modified_peptide_instance_id",
            "peptide_instance_id",
            "modified_peptide_sequence",
            "modification_mass",
            "n_observations",
        ],
        [
            {
                "modified_peptide_instance_id": "1001",
                "peptide_instance_id": "1",
                "modified_peptide_sequence": "AAAAAS[Phospho]AAAAA",
                "modification_mass": "79.9663",
                "n_observations": "50",
            },
        ],
    )

    zip_path = tmp_path / "atlas_build_canonical_test.tsv.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        for tsv_file in tables_dir.glob("*.tsv"):
            zf.write(tsv_file, tsv_file.name)

    sites = parse_peptideatlas_zip(str(zip_path))

    # Only canonical protein P04637 should have sites, not Q99999
    accessions = {s["accession"] for s in sites}
    assert "P04637" in accessions
    assert "Q99999" not in accessions
    assert len(sites) == 1
    assert sites[0]["position"] == 315


# ---------- Test 9: URL construction ----------


def test_url_no_doubled_phospho():
    """from_build() URL should not contain doubled /phospho/."""
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        PeptideAtlasPhosphoDataset,
    )

    dataset = PeptideAtlasPhosphoDataset.from_build("202512", "606")
    assert "/phospho/phospho/" not in dataset.zip_url
    assert "202512/atlas_build_606.tsv.zip" in dataset.zip_url


@pytest.mark.parametrize(
    ("table", "replacement"),
    [
        (
            "modified_peptide_instance.tsv",
            ("modified_peptide_sequence", "renamed_sequence"),
        ),
        (
            "modified_peptide_instance.tsv",
            ("n_observations", "renamed_observations"),
        ),
    ],
)
def test_parse_rejects_missing_required_columns(
    mock_pa_zip, tmp_path, table, replacement
):
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        parse_peptideatlas_zip,
    )

    path = _rewrite_zip_table(mock_pa_zip, tmp_path, table, [replacement])
    with pytest.raises(ValueError, match="missing required column"):
        parse_peptideatlas_zip(str(path))


def test_parse_rejects_empty_canonical_table(tmp_path):
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        parse_peptideatlas_zip,
    )

    original = tmp_path / "canonical.zip"
    tables = tmp_path / "canonical-tables"
    tables.mkdir()
    _write_tsv(
        tables / "protein_identification.tsv",
        ["biosequence_id", "presence_level_id"],
        [{"biosequence_id": "1", "presence_level_id": "3"}],
    )
    # The remaining required tables are valid, but have no data rows.
    _write_tsv(
        tables / "biosequence.tsv",
        [
            "biosequence_id",
            "biosequence_accession",
            "biosequence_gene_name",
            "biosequence_seq",
        ],
        [],
    )
    _write_tsv(
        tables / "peptide_instance.tsv", ["peptide_instance_id", "n_observations"], []
    )
    _write_tsv(
        tables / "peptide_mapping.tsv",
        ["peptide_instance_id", "matched_biosequence_id", "start_in_biosequence"],
        [],
    )
    _write_tsv(
        tables / "modified_peptide_instance.tsv",
        ["peptide_instance_id", "modified_peptide_sequence", "n_observations"],
        [],
    )
    with zipfile.ZipFile(original, "w") as zf:
        for path in tables.glob("*.tsv"):
            zf.write(path, path.name)

    with pytest.raises(ValueError, match="no canonical proteins"):
        parse_peptideatlas_zip(str(original))


def test_parse_normalizes_backslash_n_null_marker(mock_pa_zip, tmp_path):
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        parse_peptideatlas_zip,
    )

    path = _rewrite_zip_table(
        mock_pa_zip, tmp_path, "biosequence.tsv", [("TP53", r"\N")]
    )
    sites = parse_peptideatlas_zip(str(path))
    assert all(site["gene_symbol"] != r"\N" for site in sites)
    assert any(site["gene_symbol"] == "" for site in sites)


def test_parse_raises_when_no_phospho_sites(mock_pa_zip, tmp_path):
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        parse_peptideatlas_zip,
    )

    path = _rewrite_zip_table(
        mock_pa_zip,
        tmp_path,
        "modified_peptide_instance.tsv",
        [("[Phospho]", "[Oxidation]"), ("[167]", "[Oxidation]")],
    )
    with pytest.raises(ValueError, match="No phospho sites found"):
        parse_peptideatlas_zip(str(path))


# ---------- Test 10: _find_table_in_zip substring disambiguation ----------


def test_find_table_prefers_shorter_name_over_longer_substring_match(tmp_path):
    """A lookup for ``peptide_instance.tsv`` must not return
    ``atlas_build_606_modified_peptide_instance.tsv`` just because it also
    contains ``peptide_instance.tsv`` as a substring. With the modified table
    listed first, the old plain substring scan returned it instead of the
    real ``atlas_build_606_peptide_instance.tsv`` (#427)."""
    from hvantk.skills.peptideatlas.phospho.shared.datasets import (
        _find_table_in_zip,
    )

    zip_path = tmp_path / "prefixed.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("atlas_build_606_modified_peptide_instance.tsv", "")
        zf.writestr("atlas_build_606_peptide_instance.tsv", "")

    with zipfile.ZipFile(zip_path) as zf:
        found = _find_table_in_zip(zf, "peptide_instance.tsv")

    assert found == "atlas_build_606_peptide_instance.tsv"

    # The file name decides, not the folder depth, and the modified table is never
    # taken for a missing peptide_instance table.
    foldered = tmp_path / "foldered.zip"
    with zipfile.ZipFile(foldered, "w") as zf:
        zf.writestr("a/atlas_build_606_modified_peptide_instance.tsv", "")
        zf.writestr("long_folder/atlas_build_606_peptide_instance.tsv", "")
    with zipfile.ZipFile(foldered) as zf:
        assert (
            _find_table_in_zip(zf, "peptide_instance.tsv")
            == "long_folder/atlas_build_606_peptide_instance.tsv"
        )
    only_modified = tmp_path / "only_modified.zip"
    with zipfile.ZipFile(only_modified, "w") as zf:
        zf.writestr("atlas_build_606_modified_peptide_instance.tsv", "")
    with zipfile.ZipFile(only_modified) as zf:
        assert _find_table_in_zip(zf, "peptide_instance.tsv") is None
