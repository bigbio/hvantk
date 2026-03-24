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

from hvantk.ptm.mapper import (
    CodonMapping,
    GTFData,
    parse_ensembl_gtf,
    map_residue_to_genomic,
    map_protein_sites,
    resolve_transcript,
)
from hvantk.ptm.pipeline import (
    PTMBuildConfig,
    PTMBuildResult,
    ptm_build_pipeline,
    map_ptm_sites,
    download_ensembl_gtf,
    download_uniprot_ptm,
)
from hvantk.ptm.annotate import (
    annotate_variants_with_ptm,
)
from hvantk.ptm.analysis import (
    PTMLandscapeResult,
    PTMPopulationResult,
    ptm_landscape,
    ptm_population,
    export_ptm_strata,
)
from hvantk.ptm.plot import (
    plot_landscape_summary,
    plot_overlap_by_category,
    plot_distance_distribution,
    plot_population_af,
    encode_figure_to_base64,
)
from hvantk.ptm.report import (
    generate_report,
)
from hvantk.ptm.constants import (
    UNIPROT_API_URL,
    UNIPROT_API_FIELDS,
    UNIPROT_HUMAN_PTM_QUERY,
    ENSEMBL_GTF_URL,
    PTM_TYPE_CATEGORIES,
    DEFAULT_FLANKING_CODONS,
    PTM_OUTPUT_COLUMNS,
)

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
