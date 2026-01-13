# HGC: Hail-based Genotype Combiner

HGC is a module within hvantk that provides high-performance tools for joint genotyping workflows using [Hail](https://hail.is/). It enables efficient combination of genomic variant call format (GVCF) files, conversion between different Hail data formats (VDS, MatrixTable), and export to standard VCF format.

## Overview

The HGC module implements a complete joint genotyping pipeline with integrated quality control:

### Primary Functionality: Joint Genotyping
1. **GVCF Combination** - Combine multiple single-sample GVCF files into a unified Variant DataSet (VDS)
2. **VDS Operations** - Merge multiple VDS datasets and convert between formats
3. **MatrixTable Processing** - Convert VDS to analysis-ready MatrixTable format
4. **VCF Export** - Export processed data back to standard VCF format

### Additional Functionality: Quality Control & Visualization
5. **QC Metrics Computation** - Comprehensive sample and variant quality assessment on combined cohorts
6. **QC Visualization** - Static and interactive plots for quality control analysis
7. **QC Reports** - Professional HTML reports with embedded plots and recommendations
8. **QC-based Filtering** - Quality-based sample and variant filtering tools

## Key Features

### Core Joint Genotyping Features
- **Scalable Joint Genotyping**: Efficiently combine thousands of GVCF files using Hail's optimized combiner
- **Flexible Data Formats**: Work with VDS (storage-optimized) and MatrixTable (analysis-optimized) formats
- **Multi-allelic Handling**: Automatic splitting and normalization of multi-allelic variants
- **Genotype Annotation**: Built-in support for adjusted genotype annotations

### Quality Control Features (Post-Combination)
- **Comprehensive QC Metrics**: Sample and variant-level quality assessment for combined cohorts
- **Interactive Visualizations**: Static matplotlib and interactive Plotly plots for data exploration
- **Professional Reports**: HTML reports with embedded plots, metrics, and recommendations
- **Quality-based Filtering**: Threshold-based sample and variant filtering tools

### Interface Options
- **CLI and Python API**: Use via command-line interface or directly in Python scripts

## Installation

HGC is part of the hvantk package. Install using Poetry:

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
poetry install
poetry shell
```

## Quick Start

### Command-Line Interface

#### Core Joint Genotyping Commands

The HGC module provides four main genotype combination commands accessible via `hvantk hgc`:

```bash
# View available commands
hvantk hgc --help

# Combine GVCF files
hvantk hgc gvcf-combine -g /path/to/gvcfs -o combined.vds

# Combine VDS datasets
hvantk hgc vds-combine -i /path/to/vds_dir -o merged.vds

# Convert VDS to MatrixTable
hvantk hgc vds2mt -i combined.vds -o analysis.mt

# Export MatrixTable to VCF
hvantk hgc mt2vcf -i analysis.mt -o results.vcf.gz
```

#### Quality Control Commands (Post-Combination)

Additional QC commands for analyzing combined cohorts:

```bash
# Compute QC metrics for combined cohort
hvantk hgc compute-qc -i analysis.mt -o analysis_qc.mt

# Generate QC visualizations
hvantk hgc plot-qc -i analysis_qc.mt -o plots/ --plot-type dashboard

# Create interactive QC plots
hvantk hgc plot-qc -i analysis_qc.mt -o plots/ --interactive

# Generate comprehensive QC report
hvantk hgc qc-report -i analysis_qc.mt -o qc_report.html

# Filter based on QC metrics
hvantk hgc filter-qc -i analysis_qc.mt -o filtered.mt --min-sample-call-rate 0.95
```

### Python API

#### Core Joint Genotyping Functions

Use HGC functions directly in Python:

```python
from hvantk.hgc import (
    combine_gvcfs,
    combine_vdses,
    convert_vds_to_mt,
    convert_mt_to_multi_sample_vcf
)

# Combine GVCF files into VDS
combine_gvcfs(
    gvcf_dir="/path/to/gvcfs",
    vds_output_path="output.vds",
    tmp_path="/tmp/hail",
    save_path="combiner_plan.json",
    vdses=[],
    kwargs={}
)

# Convert VDS to MatrixTable
convert_vds_to_mt(
    vds_path="output.vds",
    output_path="analysis.mt",
    adjust_genotypes=True,
    skip_split_multi=False,
    convert_lgt_to_gt=True
)
```

#### Quality Control Functions (Post-Combination)

Additional QC functionality for combined cohorts:

```python
from hvantk.hgc import compute_full_qc, filter_samples_by_qc, filter_variants_by_qc

# Load combined MatrixTable
import hail as hl
mt = hl.read_matrix_table("analysis.mt")

# Compute comprehensive QC metrics
qc_results = compute_full_qc(mt)

# Generate visualizations
qc_results.plot_dashboard(save_path='qc_dashboard.png')
qc_results.plot_interactive_dashboard().show()
qc_results.generate_html_report('qc_report.html')

# Apply quality filters
mt_filtered = filter_samples_by_qc(
    qc_results.mt,
    min_call_rate=0.95,
    min_ti_tv_ratio=1.8
)

mt_filtered = filter_variants_by_qc(
    mt_filtered,
    min_call_rate=0.90,
    min_hwe_pvalue=1e-6
)
```

## Core Components

### 1. Combiners (`hvantk.hgc.combiners`)

Functions for combining genomic datasets:

- **`combine_gvcfs()`** - Combine GVCF files and/or existing VDS datasets into a new VDS
- **`combine_vdses()`** - Merge multiple VDS directories into a single VDS

**Advanced functions** (require direct import from `hvantk.hgc.combiners`):
- **`combine_matrix_table_rows()`** - Combine MatrixTables by rows (variants)
- **`combine_matrix_table_cols()`** - Combine MatrixTables by columns (samples)

### 2. Converters (`hvantk.hgc.converters`)

Functions for format conversion:

- **`convert_vds_to_mt()`** - Convert VDS to dense MatrixTable format
- **`convert_mt_to_multi_sample_vcf()`** - Export MatrixTable to multi-sample VCF

### 3. File Utilities (`hvantk.hgc.file_utils`)

Helper functions for file handling:

- **`validate_vcfs_paths()`** - Validate GVCF files and their indices
- **`validate_vds_paths()`** - Validate VDS directory paths
- **`check_path_exists_and_readable()`** - Verify file accessibility
- **`compress_files()`** / **`decompress_files()`** - Handle file compression
- **`sort_mts_cols()`** - Sort MatrixTable columns to match reference order

## Detailed Usage

### GVCF Combination

Combine single-sample GVCF files into a joint-called VDS:

**CLI:**
```bash
hvantk hgc gvcf-combine \
  --gvcf-dir /data/gvcfs \
  --output cohort.vds \
  --temp-dir /tmp/hail \
  --save-path combiner_plan.json
```

**Python:**
```python
from hvantk.hgc import combine_gvcfs

combine_gvcfs(
    gvcf_dir="/data/gvcfs",
    vds_output_path="cohort.vds",
    tmp_path="/tmp/hail",
    save_path="combiner_plan.json",
    vdses=[],  # Optional: existing VDS to include
    kwargs={
        # Optional interval parameters (mutually exclusive):
        # 'use_genome_default_intervals': True,
        # 'use_exome_default_intervals': True,
        # 'intervals': [list of intervals]
    },
    reference_genome="GRCh38"
)
```

**Parameters:**
- `gvcf_dir`: Directory containing `.g.vcf.gz` files with `.tbi` indices
- `vds_output_path`: Output path for the combined VDS
- `tmp_path`: Temporary directory for intermediate files
- `save_path`: Path to save the combiner execution plan (JSON)
- `vdses`: List of existing VDS paths to combine with GVCFs
- `kwargs`: Additional parameters for `hail.vds.new_combiner()`
- `reference_genome`: Reference genome build (default: "GRCh38")

**Requirements:**
- GVCF files must have corresponding `.tbi` index files
- Sufficient temporary storage (typically 2-3x input size)
- Adequate memory for Hail operations

### VDS Combination

Merge multiple VDS datasets:

**CLI:**
```bash
hvantk hgc vds-combine \
  --input-dir /data/vds_datasets \
  --output merged.vds \
  --validate \
  --overwrite
```

**Python:**
```python
from hvantk.hgc import combine_vdses

combine_vdses(
    vdses_dir="/data/vds_datasets",
    output_path="merged.vds",
    validate=True,
    overwrite=False
)
```

**Parameters:**
- `vdses_dir`: Directory containing VDS subdirectories (each ending in `.vds`)
- `output_path`: Output path for merged VDS
- `validate`: Whether to validate the combined VDS using Hail's validator
- `overwrite`: Whether to overwrite existing output

### VDS to MatrixTable Conversion

Convert storage-optimized VDS to analysis-ready MatrixTable:

**CLI:**
```bash
hvantk hgc vds2mt \
  --input cohort.vds \
  --output analysis.mt \
  --adjust-genotypes \
  --convert-lgt-to-gt \
  --overwrite
```

**Python:**
```python
from hvantk.hgc import convert_vds_to_mt

convert_vds_to_mt(
    vds_path="cohort.vds",
    output_path="analysis.mt",
    adjust_genotypes=True,
    skip_split_multi=False,
    convert_lgt_to_gt=True,
    skip_keying_by_cols=False,
    overwrite=False
)
```

**Parameters:**
- `vds_path`: Input VDS path
- `output_path`: Output MatrixTable path
- `adjust_genotypes`: Annotate with adjusted genotypes using gnomAD quality filters (requires `gnomad` package)
- `skip_split_multi`: Skip splitting multi-allelic variants (not recommended)
- `convert_lgt_to_gt`: Convert LGT (local genotype) to GT (global genotype) after splitting (recommended)
- `skip_keying_by_cols`: Skip keying MatrixTable by sample column
- `overwrite`: Whether to overwrite existing output

**Important Notes:**
- LGT to GT conversion must happen **after** splitting multi-allelic variants
- Adjusted genotype annotation requires the `gnomad` package: `pip install gnomad`
- Setting `skip_split_multi=True` and `convert_lgt_to_gt=True` will raise an error

### MatrixTable to VCF Export

Export processed MatrixTable to standard VCF format:

**CLI:**
```bash
hvantk hgc mt2vcf \
  --input analysis.mt \
  --output results.vcf.gz \
  --filter-adj \
  --min-ac 1 \
  --split-multi
```

**Python:**
```python
from hvantk.hgc import convert_mt_to_multi_sample_vcf

convert_mt_to_multi_sample_vcf(
    mt_path="analysis.mt",
    vcf_path="results.vcf.gz",
    filter_adj_genotypes=True,
    min_ac=1,
    split_multi=True
)
```

**Parameters:**
- `mt_path`: Input MatrixTable path
- `vcf_path`: Output VCF file path
- `filter_adj_genotypes`: Filter to adjusted genotypes only (recommended)
- `min_ac`: Minimum alternate allele count for variant inclusion
- `split_multi`: Split multi-allelic variants (if not already split)

**VCF Info Fields:**
The output VCF includes standard INFO fields:
- `AF`: Allele frequency
- `AC`: Allele count
- `AN`: Total number of alleles
- `call_rate`: Variant call rate

## Typical Workflow

A complete joint genotyping and QC workflow using HGC:

### Core Joint Genotyping Pipeline

```bash
# Step 1: Combine individual GVCF files
hvantk hgc gvcf-combine \
  -g /data/gvcfs \
  -o cohort.vds \
  --temp-dir /tmp/hail

# Step 2: Convert to MatrixTable for analysis
hvantk hgc vds2mt \
  -i cohort.vds \
  -o cohort.mt \
  --adjust-genotypes

# Step 3: Export results to VCF
hvantk hgc mt2vcf \
  -i cohort.mt \
  -o cohort_joint_called.vcf.gz \
  --filter-adj \
  --min-ac 2
```

### Optional: Post-Combination Quality Control

```bash
# Step 4: Compute QC metrics for combined cohort
hvantk hgc compute-qc \
  -i cohort.mt \
  -o cohort_qc.mt

# Step 5: Generate QC report
hvantk hgc qc-report \
  -i cohort_qc.mt \
  -o cohort_qc_report.html

# Step 6: Filter based on QC (optional)
hvantk hgc filter-qc \
  -i cohort_qc.mt \
  -o cohort_filtered.mt \
  --min-sample-call-rate 0.95 \
  --min-variant-call-rate 0.90

# Step 7: Export filtered results
hvantk hgc mt2vcf \
  -i cohort_filtered.mt \
  -o cohort_filtered.vcf.gz \
  --filter-adj
```

## Data Formats

### GVCF (Genomic VCF)
- Single-sample variant calls with reference blocks
- Format: `.g.vcf.gz` with `.tbi` index
- Input format for joint genotyping

### VDS (Variant DataSet)
- Hail's storage-optimized format for genomic data
- Separates reference data from variant data
- Efficient for large cohorts
- Directory structure ending in `.vds`

### MatrixTable (MT)
- Hail's analysis-optimized format
- Dense matrix representation
- Rows: variants, Columns: samples
- Directory structure ending in `.mt`

### VCF (Variant Call Format)
- Standard genomic variant format
- Compatible with most genomic tools
- Format: `.vcf.gz` with optional `.tbi` index

## Performance Considerations

### Memory Requirements
- GVCF combination: ~8-16 GB for small cohorts (<100 samples)
- Large cohorts (>1000 samples): 32-64 GB recommended
- Set Hail memory with environment variable: `export PYSPARK_SUBMIT_ARGS="--driver-memory 32g"`

### Storage Requirements
- Temporary space: 2-3x input GVCF size
- VDS output: ~50-70% of input GVCF size
- MatrixTable: ~100-150% of VDS size

### Optimization Tips
1. Use SSD storage for temporary files
2. Adjust Hail partitioning for your cluster size
3. Process intervals in parallel for large cohorts
4. Use genome/exome default intervals for better performance

## Common Options

### Interval Selection
When combining GVCFs, you can specify genomic intervals:

```python
# Use default genome intervals (recommended for WGS)
kwargs = {'use_genome_default_intervals': True}

# Use default exome intervals (recommended for WES)
kwargs = {'use_exome_default_intervals': True}

# Specify custom intervals
kwargs = {'intervals': [
    hl.parse_locus_interval('chr1:1-1000000'),
    hl.parse_locus_interval('chr2:5000000-6000000')
]}
```

### Reference Genomes
Supported reference genomes:
- `GRCh38` (default) - Human genome build 38
- `GRCh37` - Human genome build 37

```python
combine_gvcfs(..., reference_genome="GRCh38")
```

## Troubleshooting

### Common Issues

**Issue: "gnomAD package not found"**
```
Solution: Install gnomAD package or disable adjusted genotypes:
pip install gnomad
# OR
convert_vds_to_mt(..., adjust_genotypes=False)
```

**Issue: "Cannot convert LGT to GT when skip_split_multi=True"**
```
Solution: Either enable multi-allelic splitting or disable LGT conversion:
convert_vds_to_mt(..., skip_split_multi=False)
# OR
convert_vds_to_mt(..., convert_lgt_to_gt=False)
```

**Issue: "Out of memory during GVCF combination"**
```
Solution: Increase driver memory:
export PYSPARK_SUBMIT_ARGS="--driver-memory 64g pyspark-shell"
```

**Issue: "TBI index file not found"**
```
Solution: Create index files for GVCF files:
tabix -p vcf input.g.vcf.gz
```

## Examples

### Example 1: Basic Joint Genotyping Pipeline

```python
import hail as hl
from hvantk.hgc import combine_gvcfs, convert_vds_to_mt, convert_mt_to_multi_sample_vcf

# Initialize Hail
hl.init()

# Combine GVCFs
combine_gvcfs(
    gvcf_dir="data/gvcfs",
    vds_output_path="output/cohort.vds",
    tmp_path="tmp",
    save_path="output/plan.json",
    vdses=[],
    kwargs={'use_genome_default_intervals': True}
)

# Convert to MatrixTable
convert_vds_to_mt(
    vds_path="output/cohort.vds",
    output_path="output/cohort.mt",
    adjust_genotypes=True
)

# Export to VCF
convert_mt_to_multi_sample_vcf(
    mt_path="output/cohort.mt",
    vcf_path="output/cohort.vcf.gz",
    filter_adj_genotypes=True,
    min_ac=1
)
```

### Example 2: Incremental Cohort Building

```python
from hvantk.hgc import combine_gvcfs, combine_vdses

# First batch
combine_gvcfs(
    gvcf_dir="data/batch1",
    vds_output_path="output/batch1.vds",
    tmp_path="tmp",
    save_path="output/batch1_plan.json",
    vdses=[],
    kwargs={}
)

# Second batch
combine_gvcfs(
    gvcf_dir="data/batch2",
    vds_output_path="output/batch2.vds",
    tmp_path="tmp",
    save_path="output/batch2_plan.json",
    vdses=[],
    kwargs={}
)

# Merge batches
combine_vdses(
    vdses_dir="output",  # Contains batch1.vds and batch2.vds
    output_path="output/merged.vds",
    validate=True,
    overwrite=False
)
```

### Example 3: Custom Quality Filtering

```python
import hail as hl
from hvantk.hgc import convert_vds_to_mt

# Convert VDS to MT
convert_vds_to_mt(
    vds_path="cohort.vds",
    output_path="cohort.mt",
    adjust_genotypes=True
)

# Read MT and apply custom filters
mt = hl.read_matrix_table("cohort.mt")

# Filter to high-quality variants
mt = mt.filter_rows(
    (mt.info.AC > 2) &  # At least 2 alternate alleles
    (mt.variant_qc.call_rate > 0.95)  # 95% call rate
)

# Filter to high-quality samples
mt = mt.filter_cols(
    mt.sample_qc.call_rate > 0.90  # 90% sample call rate
)

# Write filtered MT
mt.write("cohort_filtered.mt")
```

## API Reference

### Main Functions

#### `combine_gvcfs(gvcf_dir, vds_output_path, tmp_path, save_path, vdses, kwargs, reference_genome='GRCh38')`
Combine GVCF files into a VDS using Hail's GVCF combiner.

**Parameters:**
- `gvcf_dir` (str): Directory containing GVCF files
- `vds_output_path` (str): Output VDS path
- `tmp_path` (str): Temporary directory path
- `save_path` (str): Path to save combiner plan
- `vdses` (List[str]): List of existing VDS paths to combine
- `kwargs` (Dict): Additional parameters for Hail's combiner
- `reference_genome` (str): Reference genome (default: 'GRCh38')

**Returns:** None

#### `combine_vdses(vdses_dir, output_path, validate=True, overwrite=False)`
Combine multiple VDS directories into a single VDS.

**Parameters:**
- `vdses_dir` (str): Directory containing VDS subdirectories
- `output_path` (str): Output merged VDS path
- `validate` (bool): Validate combined VDS (default: True)
- `overwrite` (bool): Overwrite existing output (default: False)

**Returns:** None

#### `convert_vds_to_mt(vds_path, output_path, adjust_genotypes=True, skip_split_multi=False, convert_lgt_to_gt=True, skip_keying_by_cols=False, overwrite=False)`
Convert a VDS to dense MatrixTable format.

**Parameters:**
- `vds_path` (str): Input VDS path
- `output_path` (str): Output MatrixTable path
- `adjust_genotypes` (bool): Annotate with adjusted genotypes (default: True)
- `skip_split_multi` (bool): Skip splitting multi-allelic variants (default: False)
- `convert_lgt_to_gt` (bool): Convert LGT to GT after splitting (default: True)
- `skip_keying_by_cols` (bool): Skip column keying (default: False)
- `overwrite` (bool): Overwrite existing output (default: False)

**Returns:** None

#### `convert_mt_to_multi_sample_vcf(mt_path, vcf_path, filter_adj_genotypes=True, min_ac=1, split_multi=True)`
Convert a MatrixTable to multi-sample VCF format.

**Parameters:**
- `mt_path` (str): Input MatrixTable path
- `vcf_path` (str): Output VCF file path
- `filter_adj_genotypes` (bool): Filter to adjusted genotypes (default: True)
- `min_ac` (int): Minimum alternate allele count (default: 1)
- `split_multi` (bool): Split multi-allelic variants (default: True)

**Returns:** None

### Utility Functions

#### `validate_vcfs_paths(directory, pattern=None)`
Retrieve and validate GVCF file paths in a directory.

**Parameters:**
- `directory` (str): Directory to search for GVCF files
- `pattern` (str): Glob pattern for matching files (default: None)

**Returns:** List[str] - List of validated GVCF file paths

#### `validate_vds_paths(vdses)`
Validate VDS directory paths.

**Parameters:**
- `vdses` (Union[str, List[str]]): Directory or list of VDS paths

**Returns:** List[str] - List of validated VDS paths

#### `check_path_exists_and_readable(path)`
Check if a file or directory exists and is readable.

**Parameters:**
- `path` (str): Path to check

**Returns:** str - The validated path

**Raises:** FileNotFoundError, PermissionError

#### `sort_mts_cols(mts, ref_index=0)`
Sort the column order of MatrixTables to match a reference.

**Parameters:**
- `mts` (List[hl.MatrixTable]): List of MatrixTables to sort
- `ref_index` (int): Index of reference MatrixTable (default: 0)

**Returns:** List[hl.MatrixTable] - Sorted MatrixTables

### Advanced Functions

These functions require direct import from `hvantk.hgc.combiners`:

#### `combine_matrix_table_rows(mt_paths, output_path, n_partitions, force_sort_cols=False, overwrite=False, kwargs=None)`
Combine multiple MatrixTables by rows (variants).

**Parameters:**
- `mt_paths` (List[str]): List of MatrixTable paths
- `output_path` (str): Output path for combined MatrixTable
- `n_partitions` (int): Number of partitions for output
- `force_sort_cols` (bool): Sort columns before combining (default: False)
- `overwrite` (bool): Overwrite existing output (default: False)
- `kwargs` (dict): Additional arguments for Hail's union_rows

**Returns:** None

#### `combine_matrix_table_cols(mt_paths, output_path, n_partitions, overwrite=False, kwargs=None)`
Combine multiple MatrixTables by columns (samples).

**Parameters:**
- `mt_paths` (List[str]): List of MatrixTable paths
- `output_path` (str): Output path for combined MatrixTable
- `n_partitions` (int): Number of partitions for output
- `overwrite` (bool): Overwrite existing output (default: False)
- `kwargs` (dict): Additional arguments for Hail's union_cols

**Returns:** None

## Constants

The module defines commonly used constants in `hvantk.hgc.constants`:

**Reference Genomes:**
- `HG38_GENOME_REFERENCE` - "GRCh38"
- `HG37_GENOME_REFERENCE` - "GRCh37"

**File Extensions:**
- `GVCF_EXTENSION` - ".g.vcf.gz"
- `GVCF_EXTENSION_TBI` - ".g.vcf.gz.tbi"
- `VCF_EXTENSION` - ".vcf.gz"
- `VCF_EXTENSION_TBI` - ".vcf.gz.tbi"
- `VDS_EXTENSION` - ".vds"

**Entry Fields:**
- `GT_FIELD` - "GT" (Genotype)
- `AD_FIELD` - "AD" (Allelic depths)
- `DP_FIELD` - "DP" (Read depth)
- `GQ_FIELD` - "GQ" (Genotype quality)
- `PL_FIELD` - "PL" (Phred-scaled likelihoods)
- `PID_FIELD` - "PID" (Physical phasing ID)
- `SB_FIELD` - "SB" (Strand bias)
- `MIN_DP_FIELD` - "MIN_DP" (Minimum read depth)
- `ADJ_GT_FIELD` - "adj" (Adjusted genotype flag)

**Note:** MatrixTable files use the `.mt` extension, but this is a directory structure convention rather than a defined constant in the module.

## Testing

Run HGC tests:

```bash
# Run all HGC tests
pytest hvantk/tests/hgc/ -v

# Run specific test
pytest hvantk/tests/hgc/test_gvcf_combiner.py -v
```

## Dependencies

- **hail** - Core Hail library for genomic data processing
- **gnomad** - Optional, required for adjusted genotype annotations
- **click** - For CLI interface
- **Python** ≥ 3.10

## Contributing

Contributions are welcome! Please ensure:
1. Code follows existing style patterns
2. Tests are added for new functionality
3. Documentation is updated
4. All tests pass

## References

- [Hail Documentation](https://hail.is/docs/0.2/)
- [GVCF Format Specification](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format)
- [VDS Format Details](https://hail.is/docs/0.2/vds/index.html)
- [Joint Genotyping Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)

## License

HGC is part of hvantk, released under the MIT License. See [LICENSE](../../LICENSE) for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/bigbio/hvantk/issues)
- **Documentation**: [hvantk Documentation](https://github.com/bigbio/hvantk/tree/main/docs)
- **Examples**: [hvantk Examples](https://github.com/bigbio/hvantk/tree/main/examples)
