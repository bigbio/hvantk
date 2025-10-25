"""
HGC (Hail-based Genotype Combiner) module.

This module provides functionality for combining and converting genomic variant data
using Hail, including GVCF combination, VDS operations, and MatrixTable conversions.
"""

from hvantk.hgc.combiners import combine_gvcfs, combine_vdses
from hvantk.hgc.converters import convert_vds_to_mt, convert_mt_to_multi_sample_vcf
from hvantk.hgc.file_utils import (
    check_path_exists_and_readable,
    validate_vcfs_paths,
    validate_vds_paths,
    compress_files,
    decompress_files,
    sort_mts_cols
)

__all__ = [
    # Functions from combiners
    'combine_gvcfs',
    'combine_vdses',
    # Functions from converters
    'convert_vds_to_mt',
    'convert_mt_to_multi_sample_vcf',
    # Functions from file_utils
    'check_path_exists_and_readable',
    'validate_vcfs_paths',
    'validate_vds_paths',
    'compress_files',
    'decompress_files',
    'sort_mts_cols'
]
