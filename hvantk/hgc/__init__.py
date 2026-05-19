"""
HGC (Hail-based Genotype Combiner) module.
This module provides functionality for combining and converting genomic variant data
using Hail, including GVCF combination, VDS operations, MatrixTable conversions,
comprehensive quality control analysis, and end-to-end pipeline orchestration.
"""

from hvantk.hgc.combiners import combine_gvcfs, combine_vdses
from hvantk.hgc.converters import convert_vds_to_mt, convert_mt_to_multi_sample_vcf
from hvantk.hgc.file_utils import (
    check_path_exists_and_readable,
    validate_vcfs_paths,
    validate_vds_paths,
    sort_mts_cols,
)
from hvantk.core.utils.file_utils import compress_files, decompress_files
from hvantk.hgc.qc import (
    QCMetrics,
    compute_sample_qc,
    compute_variant_qc,
    compute_full_qc,
    extract_qc_metrics,
    filter_samples_by_qc,
    filter_variants_by_qc,
    get_qc_summary_stats,
    prepare_qc_for_visualization,
    save_qc_metrics,
)
from hvantk.hgc.pipeline import (
    PipelineConfig,
    PipelineState,
    PipelineStage,
    PipelineRunner,
)

__all__ = [
    # Functions from combiners
    "combine_gvcfs",
    "combine_vdses",
    # Functions from converters
    "convert_vds_to_mt",
    "convert_mt_to_multi_sample_vcf",
    # Functions from file_utils
    "check_path_exists_and_readable",
    "validate_vcfs_paths",
    "validate_vds_paths",
    "compress_files",
    "decompress_files",
    "sort_mts_cols",
    # Classes and functions from qc
    "QCMetrics",
    "compute_sample_qc",
    "compute_variant_qc",
    "compute_full_qc",
    "extract_qc_metrics",
    "filter_samples_by_qc",
    "filter_variants_by_qc",
    "get_qc_summary_stats",
    "prepare_qc_for_visualization",
    "save_qc_metrics",
    # Classes from pipeline
    "PipelineConfig",
    "PipelineState",
    "PipelineStage",
    "PipelineRunner",
]
