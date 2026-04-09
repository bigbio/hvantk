"""PTM (Post-Translational Modification) variant classification module.

Provides tools for mapping PTM sites to genomic coordinates, cross-referencing
with genetic variants, and analyzing PTM-variant landscape across the proteome.

Example:
    >>> from hvantk.ptm import PTMBuildConfig, ptm_build_pipeline
    >>> config = PTMBuildConfig(
    ...     output_dir="data/ptm/",
    ...     output_ht="data/ptm/ptm_sites.ht",
    ...     gtf_path="data/ref/Homo_sapiens.GRCh38.113.gtf.gz",
    ...     ptm_tsv="data/ptm/uniprot-ptm-human.tsv",
    ... )
    >>> result = ptm_build_pipeline(config)
"""

import importlib as _importlib

from hvantk.ptm.constants import (
    UNIPROT_API_URL,
    UNIPROT_API_FIELDS,
    UNIPROT_HUMAN_PTM_QUERY,
    ENSEMBL_GTF_URL,
    PTM_TYPE_CATEGORIES,
    DEFAULT_FLANKING_CODONS,
    PTM_OUTPUT_COLUMNS,
)

# Lazy imports for Hail-dependent and heavy modules (PEP 562).
# Accessing any name listed here triggers on-demand loading so that
# ``from hvantk.ptm.constants import ...`` never pulls in Hail.

_LAZY_MODULES = {
    # mapper
    "CodonMapping": ("hvantk.ptm.mapper", "CodonMapping"),
    "GTFData": ("hvantk.ptm.mapper", "GTFData"),
    "parse_ensembl_gtf": ("hvantk.ptm.mapper", "parse_ensembl_gtf"),
    "map_residue_to_genomic": ("hvantk.ptm.mapper", "map_residue_to_genomic"),
    "map_protein_sites": ("hvantk.ptm.mapper", "map_protein_sites"),
    "resolve_transcript": ("hvantk.ptm.mapper", "resolve_transcript"),
    # pipeline
    "PTMBuildConfig": ("hvantk.ptm.pipeline", "PTMBuildConfig"),
    "PTMBuildResult": ("hvantk.ptm.pipeline", "PTMBuildResult"),
    "ptm_build_pipeline": ("hvantk.ptm.pipeline", "ptm_build_pipeline"),
    "map_ptm_sites": ("hvantk.ptm.pipeline", "map_ptm_sites"),
    "download_ensembl_gtf": ("hvantk.ptm.pipeline", "download_ensembl_gtf"),
    "download_uniprot_ptm": ("hvantk.ptm.pipeline", "download_uniprot_ptm"),
    # annotate (requires Hail)
    "annotate_variants_with_ptm": ("hvantk.ptm.annotate", "annotate_variants_with_ptm"),
    # analysis (requires Hail)
    "PTMLandscapeResult": ("hvantk.ptm.analysis", "PTMLandscapeResult"),
    "PTMPopulationResult": ("hvantk.ptm.analysis", "PTMPopulationResult"),
    "ptm_landscape": ("hvantk.ptm.analysis", "ptm_landscape"),
    "ptm_population": ("hvantk.ptm.analysis", "ptm_population"),
    "export_ptm_strata": ("hvantk.ptm.analysis", "export_ptm_strata"),
    # plot
    "plot_landscape_summary": ("hvantk.ptm.plot", "plot_landscape_summary"),
    "plot_overlap_by_category": ("hvantk.ptm.plot", "plot_overlap_by_category"),
    "plot_distance_distribution": ("hvantk.ptm.plot", "plot_distance_distribution"),
    "plot_population_af": ("hvantk.ptm.plot", "plot_population_af"),
    "encode_figure_to_base64": ("hvantk.ptm.plot", "encode_figure_to_base64"),
    # report
    "generate_report": ("hvantk.ptm.report", "generate_report"),
}


def __getattr__(name: str):
    if name in _LAZY_MODULES:
        mod_path, attr = _LAZY_MODULES[name]
        mod = _importlib.import_module(mod_path)
        val = getattr(mod, attr)
        globals()[name] = val  # cache for subsequent access
        return val
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    # Annotation
    "annotate_variants_with_ptm",
    # Analysis (Q1, Q3)
    "PTMLandscapeResult",
    "PTMPopulationResult",
    "ptm_landscape",
    "ptm_population",
    "export_ptm_strata",
    # Pipeline API
    "PTMBuildConfig",
    "PTMBuildResult",
    "ptm_build_pipeline",
    "map_ptm_sites",
    "download_ensembl_gtf",
    "download_uniprot_ptm",
    # Data classes
    "CodonMapping",
    "GTFData",
    # Mapper functions
    "parse_ensembl_gtf",
    "map_residue_to_genomic",
    "map_protein_sites",
    "resolve_transcript",
    # Visualization
    "plot_landscape_summary",
    "plot_overlap_by_category",
    "plot_distance_distribution",
    "plot_population_af",
    "encode_figure_to_base64",
    # Report
    "generate_report",
    # Constants
    "UNIPROT_API_URL",
    "UNIPROT_API_FIELDS",
    "UNIPROT_HUMAN_PTM_QUERY",
    "ENSEMBL_GTF_URL",
    "PTM_TYPE_CATEGORIES",
    "DEFAULT_FLANKING_CODONS",
    "PTM_OUTPUT_COLUMNS",
]
