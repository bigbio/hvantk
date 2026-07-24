"""Loading a cohort's label gene set from a GeneSetCollection JSON.

Pure Python -- no Hail. Note this uses hvantk.core.utils.gene_sets.GeneSet (the plain
named set inside a collection), NOT hvantk.core.models.GeneSet (the typed artifact).
"""
import json

import pytest

from hvantk.algorithms.cohort.labels import load_label_genes
from hvantk.algorithms.cohort.spec import CohortLabels


def _write_collection(tmp_path, sets):
    payload = {
        "gene_sets": {
            name: {
                "name": name,
                "genes": sorted(genes),
                "source": "prepare-geneset",
                "n_genes": len(genes),
                "metadata": {},
            }
            for name, genes in sets.items()
        },
        "background_genes": sorted({g for genes in sets.values() for g in genes}),
        "n_gene_sets": len(sets),
        "n_background": len({g for genes in sets.values() for g in genes}),
        "source_description": "test",
        "metadata": {"created_by": "hvantk genesets prepare"},
    }
    p = tmp_path / "panel.json"
    p.write_text(json.dumps(payload, indent=2))
    return str(p)


def test_single_set_collection_needs_no_set_name(tmp_path):
    path = _write_collection(tmp_path, {"cardiac_panel": {"MYH6", "TNNT2"}})
    genes = load_label_genes(CohortLabels(gene_set=path))
    assert genes == {"MYH6", "TNNT2"}


def test_set_name_selects_one_of_several_sets(tmp_path):
    path = _write_collection(
        tmp_path, {"cardiac": {"MYH6"}, "neuro": {"SCN1A", "STXBP1"}}
    )
    genes = load_label_genes(CohortLabels(gene_set=path, set_name="neuro"))
    assert genes == {"SCN1A", "STXBP1"}


def test_multi_set_collection_without_set_name_fails_loud_and_lists_names(tmp_path):
    path = _write_collection(tmp_path, {"cardiac": {"MYH6"}, "neuro": {"SCN1A"}})
    with pytest.raises(ValueError, match="cardiac"):
        load_label_genes(CohortLabels(gene_set=path))


def test_unknown_set_name_fails_loud(tmp_path):
    path = _write_collection(tmp_path, {"cardiac": {"MYH6"}})
    with pytest.raises(ValueError, match="nope"):
        load_label_genes(CohortLabels(gene_set=path, set_name="nope"))


def test_wrong_schema_json_fails_loud_instead_of_loading_empty(tmp_path):
    """GeneSetCollection.from_dict uses .get() throughout, so an unrelated JSON
    document loads as an EMPTY collection rather than raising. Silently empty
    labels would make every downstream AUC meaningless."""
    p = tmp_path / "wrong.json"
    p.write_text(json.dumps({"totally": "unrelated"}))
    with pytest.raises(ValueError, match="no gene sets"):
        load_label_genes(CohortLabels(gene_set=str(p)))


def test_empty_selected_set_fails_loud(tmp_path):
    path = _write_collection(tmp_path, {"cardiac": set()})
    with pytest.raises(ValueError, match="is empty"):
        load_label_genes(CohortLabels(gene_set=path))
