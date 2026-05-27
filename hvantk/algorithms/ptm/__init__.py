"""PTM (Post-Translational Modification) variant classification module.

Provides tools for mapping PTM sites to genomic coordinates, cross-referencing
with genetic variants, and analyzing PTM-variant landscape across the proteome.

For the full end-to-end workflow (including UniProt download and Hail Table build),
use :mod:`hvantk.tools.ptm.pipeline`.

Example (pure mapping core):
    >>> from hvantk.algorithms.ptm import PTMBuildConfig, ptm_build_pipeline_core
    >>> config = PTMBuildConfig(
    ...     output_dir="data/ptm/",
    ...     output_ht="data/ptm/ptm_sites.ht",
    ...     gtf_path="data/ref/Homo_sapiens.GRCh38.113.gtf.gz",
    ...     ptm_tsv="data/ptm/uniprot-ptm-human.tsv",
    ... )
    >>> result = ptm_build_pipeline_core(config)
"""

import importlib as _importlib

from hvantk.core.ptm_constants import (
    UNIPROT_API_URL,
    UNIPROT_API_FIELDS,
    UNIPROT_HUMAN_PTM_QUERY,
    ENSEMBL_GTF_URL,
    PTM_TYPE_CATEGORIES,
    DEFAULT_FLANKING_CODONS,
    PTM_OUTPUT_COLUMNS,
    PROXIMAL_BP,
    DEFAULT_MAF_THRESHOLDS,
    EXPRESSION_BIN_LABELS,
    LOG_AF_EPSILON,
    LMM_MIN_N_PTM,
    LMM_MIN_N_NONPTM,
    LMM_MIN_MIXED_GENES,
    LMM_BINNED_MIN_POS_EXPR,
    LMM_BINNED_MIN_CELL_N,
)

# Lazy imports for Hail-dependent and heavy modules (PEP 562).
# Accessing any name listed here triggers on-demand loading so that
# ``from hvantk.core.ptm_constants import ...`` never pulls in Hail.

_LAZY_MODULES = {
    # mapper
    "CodonMapping": ("hvantk.algorithms.ptm.mapper", "CodonMapping"),
    "GTFData": ("hvantk.algorithms.ptm.mapper", "GTFData"),
    "parse_ensembl_gtf": ("hvantk.algorithms.ptm.mapper", "parse_ensembl_gtf"),
    "map_residue_to_genomic": ("hvantk.algorithms.ptm.mapper", "map_residue_to_genomic"),
    "map_protein_sites": ("hvantk.algorithms.ptm.mapper", "map_protein_sites"),
    "resolve_transcript": ("hvantk.algorithms.ptm.mapper", "resolve_transcript"),
    # pipeline
    "PTMBuildConfig": ("hvantk.algorithms.ptm.pipeline", "PTMBuildConfig"),
    "PTMBuildResult": ("hvantk.algorithms.ptm.pipeline", "PTMBuildResult"),
    "ptm_build_pipeline_core": ("hvantk.algorithms.ptm.pipeline", "ptm_build_pipeline_core"),
    "map_ptm_sites": ("hvantk.algorithms.ptm.pipeline", "map_ptm_sites"),
    "download_ensembl_gtf": ("hvantk.algorithms.ptm.pipeline", "download_ensembl_gtf"),
    # annotate (requires Hail)
    "annotate_variants_with_ptm": ("hvantk.algorithms.ptm.annotate", "annotate_variants_with_ptm"),
    # analysis (requires Hail)
    "PTMLandscapeResult": ("hvantk.algorithms.ptm.analysis", "PTMLandscapeResult"),
    "PTMPopulationResult": ("hvantk.algorithms.ptm.analysis", "PTMPopulationResult"),
    "ptm_landscape": ("hvantk.algorithms.ptm.analysis", "ptm_landscape"),
    "ptm_population": ("hvantk.algorithms.ptm.analysis", "ptm_population"),
    "export_ptm_strata": ("hvantk.algorithms.ptm.analysis", "export_ptm_strata"),
    # constraint (stratified AF depletion; requires Hail at runtime)
    "PTMConstraintConfig": ("hvantk.algorithms.ptm.constraint", "PTMConstraintConfig"),
    "PTMConstraintResult": ("hvantk.algorithms.ptm.constraint", "PTMConstraintResult"),
    "run_ptm_constraint": ("hvantk.algorithms.ptm.constraint", "run_ptm_constraint"),
    "load_gene_by_group_matrix": (
        "hvantk.algorithms.ptm.constraint_expression",
        "load_gene_by_group_matrix",
    ),
    # plot
    "plot_landscape_summary": ("hvantk.algorithms.ptm.plot", "plot_landscape_summary"),
    "plot_overlap_by_category": ("hvantk.algorithms.ptm.plot", "plot_overlap_by_category"),
    "plot_distance_distribution": ("hvantk.algorithms.ptm.plot", "plot_distance_distribution"),
    "plot_population_af": ("hvantk.algorithms.ptm.plot", "plot_population_af"),
    "encode_figure_to_base64": ("hvantk.algorithms.ptm.plot", "encode_figure_to_base64"),
    # report
    "generate_report": ("hvantk.algorithms.ptm.report", "generate_report"),
    # Phase-2 atlas facade
    "PTMAtlasConfig": ("hvantk.algorithms.ptm.atlas", "PTMAtlasConfig"),
    "PTMAtlasResult": ("hvantk.algorithms.ptm.atlas", "PTMAtlasResult"),
    "build_atlas": ("hvantk.algorithms.ptm.atlas", "build_atlas"),
    # Phase-2 SYMBOL-based annotation (pandas)
    "annotate_variants_by_symbol": (
        "hvantk.algorithms.ptm.annotate",
        "annotate_variants_by_symbol",
    ),
    # Phase-2 constraint tests (statsmodels)
    "LMMResult": ("hvantk.algorithms.ptm.test", "LMMResult"),
    "BinnedLMMResult": ("hvantk.algorithms.ptm.test", "BinnedLMMResult"),
    "run_lmm": ("hvantk.algorithms.ptm.test", "run_lmm"),
    "run_binned_interaction_lmm": (
        "hvantk.algorithms.ptm.test",
        "run_binned_interaction_lmm",
    ),
    # Phase-2 report writer
    "generate_phase2_report": ("hvantk.algorithms.ptm.report", "generate_phase2_report"),
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
    # Constraint (stratified AF depletion)
    "PTMConstraintConfig",
    "PTMConstraintResult",
    "run_ptm_constraint",
    "load_gene_by_group_matrix",
    # Pipeline API
    "PTMBuildConfig",
    "PTMBuildResult",
    "ptm_build_pipeline_core",
    "map_ptm_sites",
    "download_ensembl_gtf",
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
    "PROXIMAL_BP",
    "DEFAULT_MAF_THRESHOLDS",
    "EXPRESSION_BIN_LABELS",
    "LOG_AF_EPSILON",
    "LMM_MIN_N_PTM",
    "LMM_MIN_N_NONPTM",
    "LMM_MIN_MIXED_GENES",
    "LMM_BINNED_MIN_POS_EXPR",
    "LMM_BINNED_MIN_CELL_N",
    # Phase-2 atlas
    "PTMAtlasConfig",
    "PTMAtlasResult",
    "build_atlas",
    # Phase-2 SYMBOL annotation
    "annotate_variants_by_symbol",
    # Phase-2 tests
    "LMMResult",
    "BinnedLMMResult",
    "run_lmm",
    "run_binned_interaction_lmm",
    # Phase-2 report
    "generate_phase2_report",
]
