"""Tests for ancestry constants module.

These tests verify that constants have expected values and types.
No Hail dependency required.
"""

import pytest
import importlib.util
import sys
from pathlib import Path


def _import_constants_directly():
    """Import constants module directly without going through __init__.py.

    This avoids triggering Hail import which happens in the ancestry __init__.py.
    """
    # Get the path to constants.py
    constants_path = Path(__file__).parent.parent.parent / "ancestry" / "constants.py"

    spec = importlib.util.spec_from_file_location(
        "hvantk_ancestry_constants", constants_path
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


constants = _import_constants_directly()


class TestDefaultValues:
    """Tests for default parameter values."""

    def test_min_af_in_valid_range(self):
        """min_af should be between 0 and 1."""
        assert 0 <= constants.DEFAULT_MIN_AF < 1

    def test_max_af_in_valid_range(self):
        """max_af should be between 0 and 1."""
        assert 0 < constants.DEFAULT_MAX_AF <= 1

    def test_af_range_valid(self):
        """min_af should be less than max_af."""
        assert constants.DEFAULT_MIN_AF < constants.DEFAULT_MAX_AF

    def test_call_rate_in_valid_range(self):
        """call_rate should be between 0 and 1."""
        assert 0 <= constants.DEFAULT_MIN_CALL_RATE <= 1

    def test_ld_r2_in_valid_range(self):
        """LD r2 should be between 0 and 1."""
        assert 0 < constants.DEFAULT_LD_R2 < 1

    def test_ld_window_positive(self):
        """LD window should be positive."""
        assert constants.DEFAULT_LD_WINDOW > 0

    def test_n_pcs_positive(self):
        """Number of PCs should be positive."""
        assert constants.DEFAULT_N_PCS > 0

    def test_n_pcs_classify_not_greater_than_n_pcs(self):
        """PCs for classification should not exceed total PCs."""
        assert constants.DEFAULT_N_PCS_CLASSIFY <= constants.DEFAULT_N_PCS

    def test_n_estimators_positive(self):
        """Number of RF estimators should be positive."""
        assert constants.DEFAULT_N_ESTIMATORS > 0

    def test_min_prob_in_valid_range(self):
        """min_prob should be between 0 and 1."""
        assert 0 < constants.DEFAULT_MIN_PROB <= 1

    def test_min_samples_per_pop_positive(self):
        """Minimum samples per population should be positive."""
        assert constants.MIN_SAMPLES_PER_POP > 0

    def test_shared_variants_thresholds_ordered(self):
        """MIN_SHARED_VARIANTS should be less than WARN_SHARED_VARIANTS."""
        assert constants.MIN_SHARED_VARIANTS < constants.WARN_SHARED_VARIANTS


class TestColumnNames:
    """Tests for column name constants."""

    def test_column_names_are_strings(self):
        """All column names should be strings."""
        assert isinstance(constants.ANCESTRY_COL, str)
        assert isinstance(constants.PREDICTED_ANCESTRY_COL, str)
        assert isinstance(constants.ANCESTRY_PROB_COL, str)
        assert isinstance(constants.SOURCE_COL, str)
        assert isinstance(constants.KNOWN_ANCESTRY_COL, str)

    def test_column_names_not_empty(self):
        """Column names should not be empty."""
        assert len(constants.ANCESTRY_COL) > 0
        assert len(constants.PREDICTED_ANCESTRY_COL) > 0
        assert len(constants.ANCESTRY_PROB_COL) > 0
        assert len(constants.SOURCE_COL) > 0
        assert len(constants.KNOWN_ANCESTRY_COL) > 0


class TestPopulationColors:
    """Tests for population color palette."""

    def test_population_colors_is_dict(self):
        """POPULATION_COLORS should be a dictionary."""
        assert isinstance(constants.POPULATION_COLORS, dict)

    def test_population_colors_has_all_superpops(self):
        """Should include all 1000 Genomes super-populations."""
        for superpop in constants.ALL_1KG_SUPERPOPS:
            assert superpop in constants.POPULATION_COLORS

    def test_population_colors_has_all_subpops(self):
        """Should include all 1000 Genomes sub-populations."""
        for subpop in constants.ALL_1KG_SUBPOPS:
            assert subpop in constants.POPULATION_COLORS

    def test_population_colors_has_unassigned(self):
        """Should include 'unassigned' color."""
        assert "unassigned" in constants.POPULATION_COLORS

    def test_colors_are_valid_hex(self):
        """All colors should be valid hex codes."""
        for pop, color in constants.POPULATION_COLORS.items():
            assert isinstance(color, str)
            assert color.startswith("#")
            assert len(color) == 7  # #RRGGBB format

    def test_superpop_colors_subset_of_population_colors(self):
        """SUPERPOP_COLORS should be a subset of POPULATION_COLORS."""
        for superpop, color in constants.SUPERPOP_COLORS.items():
            assert superpop in constants.POPULATION_COLORS


class TestPopulationMappings:
    """Tests for population mapping constants."""

    def test_superpop_to_subpops_has_five_superpops(self):
        """Should have 5 super-populations in mapping."""
        assert len(constants.SUPERPOP_TO_SUBPOPS) == 5

    def test_superpop_to_subpops_keys_match_all_superpops(self):
        """Keys should match ALL_1KG_SUPERPOPS."""
        assert set(constants.SUPERPOP_TO_SUBPOPS.keys()) == set(constants.ALL_1KG_SUPERPOPS)

    def test_all_subpops_count(self):
        """Should have 26 sub-populations in 1KG Phase 3."""
        assert len(constants.ALL_1KG_SUBPOPS) == 26

    def test_subpop_to_superpop_reverse_mapping(self):
        """SUBPOP_TO_SUPERPOP should be correct reverse mapping."""
        for superpop, subpops in constants.SUPERPOP_TO_SUBPOPS.items():
            for subpop in subpops:
                assert constants.SUBPOP_TO_SUPERPOP[subpop] == superpop

    def test_all_subpops_in_reverse_mapping(self):
        """All sub-populations should be in reverse mapping."""
        for subpop in constants.ALL_1KG_SUBPOPS:
            assert subpop in constants.SUBPOP_TO_SUPERPOP

    def test_afr_has_seven_subpops(self):
        """AFR should have 7 sub-populations."""
        assert len(constants.SUPERPOP_TO_SUBPOPS["AFR"]) == 7

    def test_amr_has_four_subpops(self):
        """AMR should have 4 sub-populations."""
        assert len(constants.SUPERPOP_TO_SUBPOPS["AMR"]) == 4

    def test_eas_has_five_subpops(self):
        """EAS should have 5 sub-populations."""
        assert len(constants.SUPERPOP_TO_SUBPOPS["EAS"]) == 5

    def test_eur_has_five_subpops(self):
        """EUR should have 5 sub-populations."""
        assert len(constants.SUPERPOP_TO_SUBPOPS["EUR"]) == 5

    def test_sas_has_five_subpops(self):
        """SAS should have 5 sub-populations."""
        assert len(constants.SUPERPOP_TO_SUBPOPS["SAS"]) == 5


class TestPopulationNames:
    """Tests for population full names."""

    def test_population_names_has_all_superpops(self):
        """Should have names for all super-populations."""
        for superpop in constants.ALL_1KG_SUPERPOPS:
            assert superpop in constants.POPULATION_NAMES

    def test_population_names_has_all_subpops(self):
        """Should have names for all sub-populations."""
        for subpop in constants.ALL_1KG_SUBPOPS:
            assert subpop in constants.POPULATION_NAMES

    def test_population_names_are_strings(self):
        """All population names should be non-empty strings."""
        for pop, name in constants.POPULATION_NAMES.items():
            assert isinstance(name, str)
            assert len(name) > 0


class TestKnownPopulations:
    """Tests for known population constants."""

    def test_all_1kg_superpops_is_tuple(self):
        """ALL_1KG_SUPERPOPS should be a tuple."""
        assert isinstance(constants.ALL_1KG_SUPERPOPS, tuple)

    def test_all_1kg_superpops_has_five_pops(self):
        """1000 Genomes has 5 super-populations."""
        assert len(constants.ALL_1KG_SUPERPOPS) == 5

    def test_known_hapmap_pops_is_tuple(self):
        """KNOWN_HAPMAP_POPS should be a tuple."""
        assert isinstance(constants.KNOWN_HAPMAP_POPS, tuple)

    def test_known_hapmap_pops_not_empty(self):
        """HapMap populations should not be empty."""
        assert len(constants.KNOWN_HAPMAP_POPS) > 0

    def test_backward_compatibility_alias(self):
        """KNOWN_1KG_SUPERPOPS should be alias for ALL_1KG_SUPERPOPS."""
        assert constants.KNOWN_1KG_SUPERPOPS == constants.ALL_1KG_SUPERPOPS
