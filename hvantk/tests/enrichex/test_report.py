from pathlib import Path

import pandas as pd

from hvantk.core.utils.gene_sets import GeneSet, GeneSetCollection
from hvantk.enrichex.report import generate_report


def _mock_enrichment_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_set_name": ["Microglia", "Astrocytes"],
            "n_overlap": [12, 8],
            "odds_ratio": [2.5, 1.7],
            "p_value": [1e-6, 5e-4],
            "p_adjusted": [5e-6, 1e-3],
            "significant": [True, True],
            "overlap_genes": ["APOE,TREM2", "GFAP,S100B"],
        }
    )


def _mock_burden_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_set_name": ["Microglia", "Astrocytes"],
            "odds_ratio": [1.6, 1.2],
            "ci_lower": [1.2, 0.9],
            "ci_upper": [2.1, 1.4],
            "p_value": [0.001, 0.08],
            "p_adjusted": [0.002, 0.12],
            "significant": [True, False],
        }
    )


def _write_tsv(df: pd.DataFrame, path: Path) -> Path:
    df.to_csv(path, sep="\t", index=False)
    return path


def _create_gene_sets(tmp_path: Path) -> Path:
    gs = GeneSet(name="Microglia", genes={"APOE", "TREM2"})
    collection = GeneSetCollection(
        gene_sets={gs.name: gs},
        background_genes={"APOE", "TREM2", "GFAP"},
        source_description="Test gene sets",
    )
    output = tmp_path / "gene_sets.json"
    collection.save(output)
    return output


def test_generate_report_with_inline_plots(tmp_path):
    overlap_path = _write_tsv(_mock_enrichment_df(), tmp_path / "overlap.tsv")
    burden_path = _write_tsv(_mock_burden_df(), tmp_path / "burden.tsv")
    genes_path = _create_gene_sets(tmp_path)
    output_path = tmp_path / "report.html"

    generate_report(
        output_path=str(output_path),
        overlap_results=str(overlap_path),
        burden_results=str(burden_path),
        gene_sets_path=str(genes_path),
        title="Test EnrichEx Report",
        analyst_name="QA",
        include_gene_lists=True,
        embed_static_plots=True,
    )

    contents = output_path.read_text()
    assert "Test EnrichEx Report" in contents
    assert "data:image/png;base64" in contents
    assert "Overlap Enrichment" in contents
    assert "Burden Testing" in contents
    assert "Gene Set Collection" in contents


def test_generate_report_external_plots(tmp_path):
    overlap_path = _write_tsv(_mock_enrichment_df(), tmp_path / "overlap.tsv")
    output_path = tmp_path / "report.html"

    generate_report(
        output_path=str(output_path),
        overlap_results=str(overlap_path),
        embed_static_plots=False,
        include_gene_lists=False,
        title="External Assets",
    )

    contents = output_path.read_text()
    assert "External Assets" in contents
    assert "enrichex_overlap.png" in contents
    assert (tmp_path / "enrichex_overlap.png").exists()
