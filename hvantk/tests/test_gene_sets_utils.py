# Tests for gene-set utilities in hvantk.core.utils.gene_sets

import logging

from hvantk.core.utils.gene_sets import load_sample_chd_gene_set

# Configure logging for tests
logging.basicConfig(level=logging.INFO)


class TestGeneSets:
    """Test gene-set helper functions"""

    def test_load_sample_gene_set(self):
        genes = load_sample_chd_gene_set()
        assert isinstance(genes, set)
        assert len(genes) > 0
        assert "GATA4" in genes
        assert "NKX2-5" in genes
