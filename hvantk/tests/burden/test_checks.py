import pytest
from hvantk.algorithms.burden.checks import (
    KEY_SPACES, CARRIER_MODES,
    check_key_space, check_carrier_mode, check_required_fields,
)


def test_key_spaces_match_cohort_contract():
    assert KEY_SPACES == ("gene_id", "hgnc_id", "symbol")


def test_check_key_space_rejects_unknown():
    with pytest.raises(ValueError, match="unknown gene key space 'ensembl'"):
        check_key_space("ensembl")


def test_check_carrier_mode_rejects_unknown():
    with pytest.raises(ValueError, match="unknown carrier mode 'dom'"):
        check_carrier_mode("dom")


def test_required_fields_all_present_ok():
    check_required_fields(
        row_fields=["locus", "alleles", "SYMBOL", "csq_group"],
        entry_fields=["GT"],
        col_fields=["s", "is_case"],
        gene_col="SYMBOL", route_col="csq_group", arm_col="is_case",
    )


def test_required_fields_missing_gene_names_it():
    with pytest.raises(ValueError, match="gene column 'SYMBOL' not found in MatrixTable row fields"):
        check_required_fields(
            row_fields=["locus", "alleles", "csq_group"],
            entry_fields=["GT"], col_fields=["s", "is_case"],
            gene_col="SYMBOL", route_col="csq_group", arm_col="is_case",
        )


def test_required_fields_missing_gt_names_it():
    with pytest.raises(ValueError, match="genotype entry field 'GT' not found"):
        check_required_fields(
            row_fields=["locus", "alleles", "SYMBOL", "csq_group"],
            entry_fields=["AD"], col_fields=["s", "is_case"],
            gene_col="SYMBOL", route_col="csq_group", arm_col="is_case",
        )


def test_required_fields_missing_arm_names_it():
    with pytest.raises(ValueError, match="case/control column 'is_case' not found in MatrixTable column fields"):
        check_required_fields(
            row_fields=["locus", "alleles", "SYMBOL", "csq_group"],
            entry_fields=["GT"], col_fields=["s"],
            gene_col="SYMBOL", route_col="csq_group", arm_col="is_case",
        )
