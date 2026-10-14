"""Streamer test for COSMIC CGC tissue-type gene sets.

``CosmicCGCGeneDiseaseTableStreamer.get_geneset_per_tissue()`` explodes the
multi-value ``tissue_type`` field (see ``builder.py``'s ``multi_value_fields``
and ``SKILL.md`` s4), so a gene tagged with more than one tissue code (e.g.
``"E, M"``) lands in each tissue's gene set rather than its own
combined-code group.
"""

from __future__ import annotations

import pytest

# Current (v103+) header, same column order as test_builder.py's
# legacy/partly-renamed-header tests.
_HEADERS = [
    "GENE_SYMBOL",
    "NAME",
    "COSMIC_GENE_ID",
    "CHROMOSOME",
    "GENOME_START",
    "GENOME_STOP",
    "CHR_BAND",
    "SOMATIC",
    "GERMLINE",
    "TUMOUR_TYPES_SOMATIC",
    "TUMOUR_TYPES_GERMLINE",
    "CANCER_SYNDROME",
    "TISSUE_TYPE",
    "MOLECULAR_GENETICS",
    "ROLE_IN_CANCER",
    "MUTATION_TYPES",
    "TRANSLOCATION_PARTNER",
    "OTHER_GERMLINE_MUT",
    "OTHER_SYNDROME",
    "TIER",
    "SYNONYMS",
]


def _row(gene_symbol: str, tissue_type: str) -> str:
    """Build one tab-separated data row, keyed by header name.

    Only the columns a test actually cares about are non-empty; everything
    else is blank, same convention as test_builder.py's hand-written headers.
    """
    values = {
        "GENE_SYMBOL": gene_symbol,
        "NAME": f"synthetic {gene_symbol}",
        "COSMIC_GENE_ID": "COSG99999",
        "CHROMOSOME": "1",
        "GENOME_START": "1000",
        "GENOME_STOP": "2000",
        "CHR_BAND": "1p36.33",
        "SOMATIC": "y",
        "GERMLINE": "n",
        "TISSUE_TYPE": tissue_type,
        "TIER": "1",
    }
    return "\t".join(values.get(header, "") for header in _HEADERS)


def _context():
    from hvantk.core.models.build_context import BuildContext

    return BuildContext(
        plugin="cosmic-cgc",
        dataset="cosmic-cgc:submissions",
        plugin_version="test",
        source_fingerprint="sha256:test",
        builder_commit=None,
    )


@pytest.mark.hail
def test_get_geneset_per_tissue_splits_multi_valued_codes(hail_session, tmp_path):
    """A gene tagged with two tissue codes lands in both tissue gene sets."""
    from hvantk.skills.cosmic_cgc.builder import build_cosmic_cgc_submissions
    from hvantk.skills.cosmic_cgc.streamers import CosmicCGCGeneDiseaseTableStreamer

    lines = [
        "\t".join(_HEADERS),
        _row("SYNTHMULTI", "E, M"),
        _row("SYNTHSINGLE", "E"),
    ]
    path = tmp_path / "multi-tissue.tsv"
    path.write_text("\n".join(lines) + "\n")

    artifact = build_cosmic_cgc_submissions(path, _context())
    table_path = str(tmp_path / "multi_tissue.ht")
    artifact.save(table_path)

    streamer = CosmicCGCGeneDiseaseTableStreamer(table_path, init_hail=False)
    result = streamer.get_geneset_per_tissue()

    assert set(result) == {"E", "M"}
    assert result["E"] == {"SYNTHMULTI", "SYNTHSINGLE"}
    assert result["M"] == {"SYNTHMULTI"}
