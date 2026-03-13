"""
Tests for gene set data structures, I/O, parsing, and validation.

Covers: GeneSet, GeneSetCollection (JSON + GMT round-trips), geneset_io
(detect_id_type, validate_gene_ids, parse_geneset_tsv, validate_with_hgnc),
and loading utilities (load_gene_sets, load_marker_genes).
"""

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from hvantk.utils.gene_sets import (
    GeneSet,
    GeneSetCollection,
    load_gene_sets,
    load_gene_sets_from_dict,
    load_marker_genes,
)
from hvantk.utils.geneset_io import (
    detect_id_type,
    parse_geneset_tsv,
    validate_gene_ids,
    validate_with_hgnc,
)

TESTDATA_GENESET = Path(__file__).parent.parent / "testdata" / "prepare_geneset"
TESTDATA_ENRICHEX = Path(__file__).parent / "testdata"


# ---------------------------------------------------------------------------
# geneset_io: detect_id_type
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "token,expected",
    [
        ("ENSG00000141510", "ensembl_gene"),
        ("ENSG00000141510.12", "ensembl_gene"),
        ("ENST00000269305", "ensembl_transcript"),
        ("ENST00000269305.8", "ensembl_transcript"),
        ("7157", "entrez"),
        ("Trp53", "mouse_symbol"),
        ("BRCA 1", "whitespace"),
        ("TP53", "symbol"),
        ("HLA-DRB1", "symbol"),
        ("TP", "symbol"),  # short symbols not classified as mouse
    ],
    ids=[
        "ensembl_gene", "ensembl_gene_versioned", "ensembl_transcript",
        "ensembl_transcript_versioned", "entrez", "mouse", "whitespace",
        "symbol", "symbol_hyphen", "short_symbol",
    ],
)
def test_detect_id_type(token, expected):
    assert detect_id_type(token) == expected


# ---------------------------------------------------------------------------
# geneset_io: validate_gene_ids
# ---------------------------------------------------------------------------


def test_validate_gene_ids_mixed():
    valid, problems = validate_gene_ids(
        ["BRCA1", "ENSG00000141510", "7157", "TP53"]
    )
    assert valid == ["BRCA1", "TP53"]
    assert "ensembl_gene" in problems
    assert "entrez" in problems


# ---------------------------------------------------------------------------
# geneset_io: parse_geneset_tsv
# ---------------------------------------------------------------------------


def test_parse_geneset_tsv_valid():
    result = parse_geneset_tsv(TESTDATA_GENESET / "valid_panels.tsv")
    assert set(result.gene_sets.keys()) == {"cardiac_panel", "epilepsy_panel", "ras_pathway"}
    assert result.gene_sets["cardiac_panel"] == ["MYH7", "TNNT2", "LMNA", "SCN5A", "TTN"]
    assert result.n_lines_parsed == 15
    assert result.n_duplicates == 0


@pytest.mark.parametrize(
    "content,error_match",
    [
        ("panel_a\tBRCA1\npanel_a\tTP53\textra\n", "expected 2 tab-separated columns"),
        ("panel_a\tENSG00000141510\npanel_a\tENSG00000012048\npanel_b\tBRCA1\n", "Ensembl gene IDs"),
        ("panel_a\tTrp53\npanel_a\tScn1a\npanel_a\tBrca1\n", "mouse gene symbols"),
        ("panel_a\t7157\npanel_a\t672\npanel_a\t1956\n", "Entrez Gene IDs"),
    ],
    ids=["wrong_columns", "ensembl", "mouse", "entrez"],
)
def test_parse_geneset_tsv_rejects_invalid(tmp_path, content, error_match):
    tsv = tmp_path / "bad.tsv"
    tsv.write_text(content)
    with pytest.raises(ValueError, match=error_match):
        parse_geneset_tsv(tsv)


def test_parse_geneset_tsv_deduplicates(tmp_path):
    tsv = tmp_path / "dupes.tsv"
    tsv.write_text("panel_a\tBRCA1\npanel_a\tTP53\npanel_a\tBRCA1\npanel_b\tEGFR\npanel_b\tEGFR\npanel_b\tKRAS\n")
    result = parse_geneset_tsv(tsv)
    assert result.gene_sets["panel_a"] == ["BRCA1", "TP53"]
    assert result.n_duplicates == 2
    assert any("Duplicate gene" in w for w in result.warnings)


# ---------------------------------------------------------------------------
# geneset_io: validate_with_hgnc
# ---------------------------------------------------------------------------


def _mock_hgnc_maps():
    canonical = {"BRCA1", "BRCA2", "TP53", "EGFR", "ERCC1"}
    alias_to_canonical = {"FANCD1": "BRCA2", "ERCC11": "ERCC1", "RNF53": "BRCA1"}
    canonical_to_aliases = {"BRCA2": ["FANCD1"], "ERCC1": ["ERCC11"], "BRCA1": ["RNF53"]}
    return canonical, alias_to_canonical, canonical_to_aliases


@patch("hvantk.utils.gene_aliases._load_hgnc_symbol_maps")
def test_validate_with_hgnc_alias_resolution(mock_load):
    mock_load.return_value = _mock_hgnc_maps()
    gene_sets = {"panel": ["FANCD1", "TP53", "ERCC11"]}
    vr = validate_with_hgnc(gene_sets, "/fake/hgnc.tsv")
    assert vr.aliases_resolved == {"FANCD1": "BRCA2", "ERCC11": "ERCC1"}
    assert vr.gene_sets["panel"] == ["BRCA2", "TP53", "ERCC1"]


@patch("hvantk.utils.gene_aliases._load_hgnc_symbol_maps")
def test_validate_with_hgnc_unrecognized(mock_load):
    mock_load.return_value = _mock_hgnc_maps()
    gene_sets = {"panel": ["BRCA1", "FAKEGENE"]}
    vr = validate_with_hgnc(gene_sets, "/fake/hgnc.tsv")
    assert vr.unrecognized == {"FAKEGENE"}
    assert "FAKEGENE" in vr.gene_sets["panel"]


# ---------------------------------------------------------------------------
# GeneSet + GeneSetCollection: round-trips
# ---------------------------------------------------------------------------


def test_geneset_dict_round_trip():
    gs = GeneSet("test", {"BRCA1", "TP53", "EGFR"}, source="manual", metadata={"x": 1})
    restored = GeneSet.from_dict(gs.to_dict())
    assert restored.name == gs.name
    assert restored.genes == gs.genes
    assert restored.source == gs.source
    assert restored.metadata == gs.metadata


def test_collection_json_round_trip(tmp_path):
    gs1 = GeneSet("set1", {"A", "B", "C"})
    gs2 = GeneSet("set2", {"D", "E", "F"})
    collection = GeneSetCollection(
        gene_sets={"set1": gs1, "set2": gs2},
        background_genes={"A", "B", "C", "D", "E", "F"},
        source_description="test",
    )
    path = tmp_path / "gene_sets.json"
    collection.save(path)
    loaded = GeneSetCollection.load(path)
    assert len(loaded) == 2
    assert loaded.get("set1").genes == gs1.genes
    assert loaded.get("set2").genes == gs2.genes
    assert loaded.background_genes == collection.background_genes


def test_collection_gmt_round_trip(tmp_path):
    gs1 = GeneSet("cancer_genes", {"BRCA1", "TP53", "EGFR"}, source="oncology")
    gs2 = GeneSet("cell_cycle", {"CDK1", "CDK2"}, source="cell biology")
    original = GeneSetCollection(
        gene_sets={"cancer_genes": gs1, "cell_cycle": gs2},
        background_genes={"BRCA1", "TP53", "EGFR", "CDK1", "CDK2"},
    )
    gmt_path = tmp_path / "round_trip.gmt"
    original.save_gmt(gmt_path)
    loaded = GeneSetCollection.load_gmt(gmt_path)
    assert len(loaded) == 2
    assert loaded.get("cancer_genes").genes == gs1.genes
    assert loaded.get("cell_cycle").genes == gs2.genes


# ---------------------------------------------------------------------------
# Loading utilities
# ---------------------------------------------------------------------------


def test_load_gene_sets_auto_detect(tmp_path):
    """JSON and GMT auto-detected by extension."""
    # JSON
    gs = GeneSet("test", {"A", "B"})
    col = GeneSetCollection(gene_sets={"test": gs}, background_genes={"A", "B", "C"})
    json_path = tmp_path / "gene_sets.json"
    col.save(json_path)
    assert load_gene_sets(json_path).get("test").genes == {"A", "B"}

    # GMT
    gmt_path = tmp_path / "gene_sets.gmt"
    gmt_path.write_text("set1\tdesc\tGENE1\tGENE2\n")
    assert load_gene_sets(gmt_path).get("set1").genes == {"GENE1", "GENE2"}

    # Unknown format
    unknown = tmp_path / "gene_sets.xyz"
    unknown.write_text("content")
    with pytest.raises(ValueError, match="Unknown gene set file format"):
        load_gene_sets(unknown)


def test_load_gene_sets_from_dict():
    collection = load_gene_sets_from_dict(
        {"set1": ["A", "B", "C"], "set2": ["D", "E", "F"]}, source="test"
    )
    assert len(collection) == 2
    assert collection.get("set1").n_genes == 3
    assert len(collection.background_genes) == 6


def test_load_marker_genes():
    test_file = TESTDATA_ENRICHEX / "marker_genes.tsv"
    if not test_file.exists():
        pytest.skip("Test data file not found")
    collection = load_marker_genes(test_file, cluster_column="cluster", gene_column="gene")
    assert len(collection) == 4
    assert "CD3D" in collection.get("T_cell").genes
