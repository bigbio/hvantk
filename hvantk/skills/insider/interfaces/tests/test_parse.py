"""Pure-Python parse of the INSIDER interface table (no Hail)."""

from __future__ import annotations

import pytest

FIXTURE = (
    "hvantk/skills/insider/interfaces/tests/testdata/raw/interfaces/"
    "H_sapiens_interfacesALL.txt"
)


@pytest.mark.parametrize(
    "field,expected",
    [
        ("[]", set()),
        ("", set()),
        ("[265]", {265}),
        ("[5,7-9]", {5, 7, 8, 9}),
        ("[29,33,36-38,41]", {29, 33, 36, 37, 38, 41}),
        ("[ 5 , 6 ]", {5, 6}),
        ("[5,junk,7]", {5, 7}),  # third-party dump: skip, never abort
    ],
)
def test_expand_ires(field, expected):
    from hvantk.skills.insider.interfaces.parse import expand_ires

    assert expand_ires(field) == expected


def test_parse_reduces_pairs_to_proteins_from_both_sides():
    """Every row feeds BOTH its proteins; residues union across interactions."""
    from hvantk.skills.insider.interfaces.parse import parse_interfaces

    df = parse_interfaces(FIXTURE).set_index("uniprot_id")
    assert sorted(df.index) == ["Q00001", "Q00002", "Q00003"]

    # Q00001: partners {Q00002, Q00003}; residues {5,7,8,9} u {5,20} -> 5 distinct
    # (residue 5 recurs at two interfaces and is counted once).
    assert df.loc["Q00001", "n_partners"] == 2
    assert df.loc["Q00001", "n_interface_residues"] == 5
    assert df.loc["Q00001", "n_partners_experimental"] == 1  # PDB only
    assert df.loc["Q00001", "n_partners_predicted"] == 1  # ECLAIR only

    # Q00002: appears as P2 (empty IRES) and as P1 ([11]).
    assert df.loc["Q00002", "n_partners"] == 2
    assert df.loc["Q00002", "n_interface_residues"] == 1

    # Q00003: only ever P2; both its partners are experimental.
    assert df.loc["Q00003", "n_partners"] == 2
    assert df.loc["Q00003", "n_interface_residues"] == 2  # {3} u {3,4}
    assert df.loc["Q00003", "n_partners_experimental"] == 2
    assert df.loc["Q00003", "n_partners_predicted"] == 0


def test_missing_column_is_rejected(tmp_path):
    from hvantk.skills.insider.interfaces.parse import parse_interfaces

    p = tmp_path / "bad.txt"
    p.write_text("P1\tP2\tSource\n" "Q1\tQ2\tECLAIR\n")
    with pytest.raises(ValueError, match="missing column"):
        parse_interfaces(str(p))
