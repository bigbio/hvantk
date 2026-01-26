# HGC Pipeline Command Implementation Plan

**Version:** 1.0  
**Date:** January 14, 2026  
**Status:** Planning  

---

## 1. Executive Summary

This document outlines the implementation plan for adding a new `pipeline` command to the HGC (Hail-based Genotype Combiner) CLI module. The pipeline command will orchestrate an end-to-end workflow for processing gVCF files through combination, quality control, and VCF export stages.

### Primary Use Case

```bash
hvantk hgc pipeline -i /path/to/gvcfs -o /path/to/output_dir
```

---

## 2. Current Architecture Analysis

### 2.1 Existing CLI Structure

The HGC CLI is implemented in `hvantk/commands/hgc_cli.py` using Click framework with the following command group structure:

```
hvantk hgc
├── gvcf-combine    # Combine GVCF files → VDS
├── vds-combine     # Combine VDS datasets
├── vds2mt          # Convert VDS → MatrixTable
├── mt2vcf          # Convert MatrixTable → VCF
├── compute-qc      # Compute sample/variant QC metrics
├── filter-qc       # Filter MT based on QC thresholds
├── qc-summary      # Generate QC summary statistics
├── plot-qc         # Generate QC visualization plots
└── qc-report       # Generate HTML QC report
```

### 2.2 Core Module Functions

Located in `hvantk/hgc/`:

| Module | Key Functions |
|--------|--------------|
| `combiners.py` | `combine_gvcfs()`, `combine_vdses()` |
| `converters.py` | `convert_vds_to_mt()`, `convert_mt_to_multi_sample_vcf()` |
| `qc.py` | `compute_sample_qc()`, `compute_variant_qc()`, `compute_full_qc()`, `filter_samples_by_qc()`, `filter_variants_by_qc()`, `save_qc_metrics()` |
| `file_utils.py` | `check_path_exists_and_readable()`, `validate_vcfs_paths()`, `validate_vds_paths()` |

### 2.3 Existing Utility Functions (in hgc_cli.py)

- `setup_logging_for_hgc()` - Configure logging
- `expand_file_patterns()` - Expand glob patterns
- `validate_input_files()` - Validate input paths
- `validate_output_path()` - Validate and create output directories
- `estimate_resource_requirements()` - Estimate memory/partitions

---

## 3. Pipeline Workflow Design

### 3.1 End-to-End Workflow Stages

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                         HGC Pipeline Workflow                                    │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                  │
│   Stage 1       Stage 2       Stage 3a      Stage 3b       Stage 4              │
│  ┌─────────┐  ┌─────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐            │
│  │ Combine │─▶│ VDS to  │─▶│  Sample  │─▶│ Variant  │─▶│  Export  │            │
│  │  gVCFs  │  │   MT    │  │    QC    │  │    QC    │  │  pVCF    │            │
│  └─────────┘  └─────────┘  └──────────┘  └──────────┘  └──────────┘            │
│                                                                                  │
│      ↓            ↓             ↓             ↓             ↓                    │
│   .vds         .mt         sample_qc.ht  variant_qc.ht  .vcf.bgz                │
│                                                                                  │
└─────────────────────────────────────────────────────────────────────────────────┘
```

### 3.2 Pipeline Stages Detail

| Stage | Description | Input | Output | Skip Flag |
|-------|-------------|-------|--------|-----------|
| 1. Combine gVCFs | Combine gVCF files into VDS | gVCF directory | `combined.vds` | `--skip-combine-gvcfs` |
| 2. VDS to MT | Convert VDS to MatrixTable | VDS | `combined.mt` | `--skip-vds-to-mt` |
| 3a. Sample QC | Compute sample-level QC | MT | `sample_qc.ht` | `--skip-compute-sample-qc` |
| 3b. Variant QC | Compute variant-level QC | MT | `variant_qc.ht` | `--skip-compute-variant-qc` |
| 4. Export pVCF | Export to multi-sample VCF | MT | `cohort.vcf.bgz` | `--skip-export-pvcf` |

---

## 4. Command Interface Specification

### 4.1 Primary Command Signature

```python
@hgc_group.command(name="pipeline")
@click.option('-i', '--input-dir', type=click.Path(exists=True), required=True,
              help='Path to directory containing input gVCF files')
@click.option('-o', '--output-dir', type=click.Path(), required=True,
              help='Path to output directory')
```

### 4.2 Stage Control Flags

```python
# Skip flags (default: False - run all stages)
@click.option('--skip-combine-gvcfs', is_flag=True, default=False,
              help='Skip combining gVCF files (use existing VDS)')
@click.option('--skip-vds-to-mt', is_flag=True, default=False,
              help='Skip VDS to MatrixTable conversion (use existing MT)')
@click.option('--skip-compute-sample-qc', is_flag=True, default=False,
              help='Skip computing sample QC metrics')
@click.option('--skip-compute-variant-qc', is_flag=True, default=False,
              help='Skip computing variant QC metrics')
@click.option('--skip-export-pvcf', is_flag=True, default=False,
              help='Skip exporting the cohort (project) VCF')
```

### 4.3 Configuration Flags

```python
# Path overrides
@click.option('--vds-path', type=click.Path(), default=None,
              help='Path to existing VDS (required if --skip-combine-gvcfs)')
@click.option('--mt-path', type=click.Path(), default=None,
              help='Path to existing MatrixTable (required if --skip-vds-to-mt)')

# Processing configuration
@click.option('--tmp-dir', type=click.Path(), default=None,
              help='Path to temporary directory for intermediate files')
@click.option('--reference-genome', type=click.Choice(['GRCh37', 'GRCh38']),
              default='GRCh38', help='Reference genome build')
@click.option('--n-partitions', type=int, default=None,
              help='Number of partitions for parallel processing')
@click.option('--overwrite', is_flag=True, default=False,
              help='Overwrite existing output files')
@click.option('--dry-run', is_flag=True, default=False,
              help='Show what would be done without executing')
```

### 4.4 QC Configuration Flags

```python
# QC thresholds (for filtering, optional)
@click.option('--min-sample-call-rate', type=float, default=0.85,
              help='Minimum sample call rate for QC filtering')
@click.option('--min-variant-call-rate', type=float, default=0.85,
              help='Minimum variant call rate for QC filtering')
@click.option('--apply-qc-filters', is_flag=True, default=False,
              help='Apply QC filters before exporting pVCF')
```

### 4.5 Output Control Flags

```python
# Output preferences
@click.option('--keep-intermediates', is_flag=True, default=True,
              help='Keep intermediate files (VDS, MT)')
@click.option('--generate-qc-report', is_flag=True, default=False,
              help='Generate HTML QC report after computing QC metrics')
@click.option('--output-prefix', type=str, default='cohort',
              help='Prefix for output files')
```

---

## 5. Implementation Steps

### Phase 1: Core Pipeline Infrastructure

#### Step 1.1: Create Pipeline State Manager

Create a new module `hvantk/hgc/pipeline.py` to manage pipeline state:

```python
# hvantk/hgc/pipeline.py

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict, Any
from enum import Enum
import json
import logging

class PipelineStage(Enum):
    COMBINE_GVCFS = "combine_gvcfs"
    VDS_TO_MT = "vds_to_mt"
    COMPUTE_SAMPLE_QC = "compute_sample_qc"
    COMPUTE_VARIANT_QC = "compute_variant_qc"
    EXPORT_PVCF = "export_pvcf"

@dataclass
class PipelineConfig:
    input_dir: str
    output_dir: str
    tmp_dir: Optional[str] = None
    reference_genome: str = "GRCh38"
    n_partitions: Optional[int] = None
    overwrite: bool = False
    output_prefix: str = "cohort"
    
    # Stage skip flags
    skip_combine_gvcfs: bool = False
    skip_vds_to_mt: bool = False
    skip_compute_sample_qc: bool = False
    skip_compute_variant_qc: bool = False
    skip_export_pvcf: bool = False
    
    # Path overrides
    vds_path: Optional[str] = None
    mt_path: Optional[str] = None
    
    # QC configuration
    min_sample_call_rate: float = 0.85
    min_variant_call_rate: float = 0.85
    apply_qc_filters: bool = False
    
    # Output options
    keep_intermediates: bool = True
    generate_qc_report: bool = False

@dataclass
class PipelineState:
    config: PipelineConfig
    current_stage: Optional[PipelineStage] = None
    completed_stages: list = field(default_factory=list)
    outputs: Dict[str, str] = field(default_factory=dict)
    errors: list = field(default_factory=list)
    
    def save(self, path: Path):
        """Save pipeline state to JSON."""
        ...
    
    @classmethod
    def load(cls, path: Path) -> 'PipelineState':
        """Load pipeline state from JSON."""
        ...
```

#### Step 1.2: Create Pipeline Runner

```python
# hvantk/hgc/pipeline.py (continued)

class PipelineRunner:
    """Orchestrates the HGC pipeline workflow."""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.state = PipelineState(config=config)
        self.logger = logging.getLogger(__name__)
        self._setup_output_paths()
    
    def _setup_output_paths(self):
        """Initialize output paths based on configuration."""
        output_dir = Path(self.config.output_dir)
        prefix = self.config.output_prefix
        
        self.paths = {
            'vds': self.config.vds_path or str(output_dir / f'{prefix}.vds'),
            'mt': self.config.mt_path or str(output_dir / f'{prefix}.mt'),
            'sample_qc': str(output_dir / 'qc' / f'{prefix}_sample_qc.ht'),
            'variant_qc': str(output_dir / 'qc' / f'{prefix}_variant_qc.ht'),
            'pvcf': str(output_dir / f'{prefix}.vcf.bgz'),
            'qc_report': str(output_dir / 'qc' / f'{prefix}_qc_report.html'),
            'state': str(output_dir / '.pipeline_state.json'),
        }
    
    def run(self) -> PipelineState:
        """Execute the pipeline."""
        ...
    
    def _run_stage(self, stage: PipelineStage) -> bool:
        """Execute a single pipeline stage."""
        ...
```

### Phase 2: CLI Command Implementation

#### Step 2.1: Add Pipeline Command to hgc_cli.py

Add the pipeline command after the existing commands in `hvantk/commands/hgc_cli.py`:

```python
@hgc_group.command(name="pipeline")
@click.option('-i', '--input-dir', type=click.Path(exists=True), required=True,
              help='Path to directory containing input gVCF files')
@click.option('-o', '--output-dir', type=click.Path(), required=True,
              help='Path to output directory')
# ... all options as specified in Section 4 ...
@click.pass_context
def pipeline(ctx, input_dir, output_dir, ...):
    """
    Run end-to-end gVCF processing pipeline.
    
    Orchestrates the complete workflow from gVCF files to cohort VCF:
    
    \b
    Stages:
      1. Combine gVCFs into VDS (Variant Dataset)
      2. Convert VDS to MatrixTable
      3. Compute sample QC metrics
      4. Compute variant QC metrics
      5. Export to project VCF (pVCF)
    
    \b
    Examples:
      # Run full pipeline
      hvantk hgc pipeline -i /path/to/gvcfs -o /path/to/output
      
      # Skip sample QC
      hvantk hgc pipeline -i /path/to/gvcfs -o /path/to/output --skip-compute-sample-qc
      
      # Start from existing VDS
      hvantk hgc pipeline -i /path/to/gvcfs -o /path/to/output \\
          --skip-combine-gvcfs --vds-path /path/to/existing.vds
      
      # Apply QC filters before export
      hvantk hgc pipeline -i /path/to/gvcfs -o /path/to/output \\
          --apply-qc-filters --min-sample-call-rate 0.9
    """
    try:
        from hvantk.hgc.pipeline import PipelineConfig, PipelineRunner
        
        # Create configuration
        config = PipelineConfig(
            input_dir=input_dir,
            output_dir=output_dir,
            # ... map all CLI options to config ...
        )
        
        # Validate configuration
        if config.skip_combine_gvcfs and not config.vds_path:
            click.echo("❌ --vds-path required when using --skip-combine-gvcfs", err=True)
            ctx.exit(1)
        
        if config.skip_vds_to_mt and not config.mt_path:
            click.echo("❌ --mt-path required when using --skip-vds-to-mt", err=True)
            ctx.exit(1)
        
        # Create and run pipeline
        runner = PipelineRunner(config)
        
        if dry_run:
            runner.show_plan()
            return
        
        state = runner.run()
        
        if state.errors:
            click.echo("❌ Pipeline completed with errors", err=True)
            ctx.exit(1)
        
        click.echo("✅ Pipeline completed successfully!")
        
    except Exception as e:
        logger.exception(f"Pipeline failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)
```

### Phase 3: Pipeline Stage Implementations

#### Step 3.1: Implement Stage Runners

Each stage will have its own method in `PipelineRunner`:

```python
def _run_combine_gvcfs(self) -> str:
    """Stage 1: Combine gVCFs into VDS."""
    from hvantk.hgc import combine_gvcfs
    
    click.echo("🔄 [1/5] Combining gVCF files into VDS...")
    
    combine_gvcfs(
        gvcf_dir=self.config.input_dir,
        vds_output_path=self.paths['vds'],
        tmp_path=self.config.tmp_dir or tempfile.mkdtemp(),
        save_path=f"{self.paths['vds']}.plan",
        vdses=[],
        kwargs={},
        reference_genome=self.config.reference_genome
    )
    
    click.echo(f"   ✓ VDS created: {self.paths['vds']}")
    return self.paths['vds']

def _run_vds_to_mt(self, vds_path: str) -> str:
    """Stage 2: Convert VDS to MatrixTable."""
    from hvantk.hgc import convert_vds_to_mt
    
    click.echo("🔄 [2/5] Converting VDS to MatrixTable...")
    
    convert_vds_to_mt(
        vds_path=vds_path,
        output_path=self.paths['mt'],
        adjust_genotypes=True,
        skip_split_multi=False,
        skip_validation=False,
        skip_keying_by_cols=False,
        overwrite=self.config.overwrite
    )
    
    click.echo(f"   ✓ MatrixTable created: {self.paths['mt']}")
    return self.paths['mt']

def _run_compute_sample_qc(self, mt_path: str) -> str:
    """Stage 3a: Compute sample QC metrics."""
    import hail as hl
    from hvantk.hgc import compute_sample_qc, save_qc_metrics
    
    click.echo("🔄 [3/5] Computing sample QC metrics...")
    
    mt = hl.read_matrix_table(mt_path)
    mt_qc = compute_sample_qc(mt)
    sample_qc_ht = mt_qc.cols().select('sample_qc')
    sample_qc_ht.write(self.paths['sample_qc'], overwrite=self.config.overwrite)
    
    click.echo(f"   ✓ Sample QC saved: {self.paths['sample_qc']}")
    return self.paths['sample_qc']

def _run_compute_variant_qc(self, mt_path: str) -> str:
    """Stage 3b: Compute variant QC metrics."""
    import hail as hl
    from hvantk.hgc import compute_variant_qc
    
    click.echo("🔄 [4/5] Computing variant QC metrics...")
    
    mt = hl.read_matrix_table(mt_path)
    mt_qc = compute_variant_qc(mt)
    variant_qc_ht = mt_qc.rows().select('variant_qc')
    variant_qc_ht.write(self.paths['variant_qc'], overwrite=self.config.overwrite)
    
    click.echo(f"   ✓ Variant QC saved: {self.paths['variant_qc']}")
    return self.paths['variant_qc']

def _run_export_pvcf(self, mt_path: str) -> str:
    """Stage 5: Export to project VCF."""
    from hvantk.hgc import convert_mt_to_multi_sample_vcf
    
    click.echo("🔄 [5/5] Exporting cohort VCF...")
    
    convert_mt_to_multi_sample_vcf(
        mt_path=mt_path,
        vcf_path=self.paths['pvcf'],
        filter_adj_genotypes=True,
        min_ac=1,
        split_multi=True
    )
    
    click.echo(f"   ✓ pVCF exported: {self.paths['pvcf']}")
    return self.paths['pvcf']
```

---

## 6. File Changes Summary

### 6.1 New Files

| File | Description |
|------|-------------|
| `hvantk/hgc/pipeline.py` | Pipeline orchestration logic, config, state management |
| `hvantk/hgc/tests/test_pipeline.py` | Unit tests for pipeline module |

### 6.2 Modified Files

| File | Changes |
|------|---------|
| `hvantk/hgc/__init__.py` | Export new pipeline classes |
| `hvantk/commands/hgc_cli.py` | Add `pipeline` command |

---

## 7. Output Directory Structure

```
output_dir/
├── cohort.vds/                    # Stage 1: Combined VDS
├── cohort.mt/                     # Stage 2: MatrixTable
├── cohort.vcf.bgz                 # Stage 5: Project VCF
├── cohort.vcf.bgz.tbi             # VCF index
├── qc/
│   ├── cohort_sample_qc.ht/       # Stage 3a: Sample QC Hail Table
│   ├── cohort_variant_qc.ht/      # Stage 3b: Variant QC Hail Table
│   ├── cohort_sample_qc.csv       # Sample QC metrics (CSV export)
│   ├── cohort_variant_qc.csv      # Variant QC metrics (CSV export)
│   └── cohort_qc_report.html      # Optional: HTML report
├── logs/
│   └── pipeline_YYYYMMDD_HHMMSS.log
└── .pipeline_state.json           # Pipeline state for resumption
```

---

## 8. Error Handling & Recovery

### 8.1 Pipeline State Persistence

The pipeline will save state after each successful stage, enabling:
- Resume from failure point
- Skip already completed stages
- Audit trail of processing

### 8.2 Validation Checks

| Check | When | Action on Failure |
|-------|------|-------------------|
| Input directory exists | Before stage 1 | Exit with error |
| gVCF files found | Before stage 1 | Exit with error |
| VDS exists | Before stage 2 (if skip-combine) | Exit with error |
| MT exists | Before stage 3-5 (if skip earlier) | Exit with error |
| Disk space | Before each stage | Warning or exit |
| Hail initialized | Before stage 1 | Initialize Hail |

---

## 9. Testing Strategy

### 9.1 Unit Tests

```python
# hvantk/hgc/tests/test_pipeline.py

def test_pipeline_config_defaults():
    """Test PipelineConfig has correct defaults."""
    ...

def test_pipeline_config_validation():
    """Test configuration validation logic."""
    ...

def test_pipeline_path_setup():
    """Test output path generation."""
    ...

def test_pipeline_skip_flags():
    """Test skip flag behavior."""
    ...

def test_pipeline_state_persistence():
    """Test state save/load."""
    ...
```

### 9.2 Integration Tests

```python
def test_pipeline_dry_run():
    """Test dry-run shows plan without execution."""
    ...

def test_pipeline_full_execution(sample_gvcfs):
    """Test full pipeline execution with sample data."""
    ...

def test_pipeline_resume_from_failure():
    """Test pipeline can resume from saved state."""
    ...
```

### 9.3 CLI Tests

```bash
# Test help
hvantk hgc pipeline --help

# Test dry-run
hvantk hgc pipeline -i ./test_gvcfs -o ./test_output --dry-run

# Test validation errors
hvantk hgc pipeline -i ./nonexistent -o ./output  # Should fail
hvantk hgc pipeline -i ./gvcfs -o ./output --skip-combine-gvcfs  # Should fail (no --vds-path)
```

---

## 10. Usage Examples

### 10.1 Full Pipeline

```bash
# Run complete end-to-end pipeline
hvantk hgc pipeline \
    -i /data/project/gvcfs \
    -o /data/project/output \
    --reference-genome GRCh38
```

### 10.2 Skip Stages

```bash
# Skip QC stages (just combine and export)
hvantk hgc pipeline \
    -i /data/project/gvcfs \
    -o /data/project/output \
    --skip-compute-sample-qc \
    --skip-compute-variant-qc
```

### 10.3 Resume from VDS

```bash
# Start from existing VDS
hvantk hgc pipeline \
    -i /data/project/gvcfs \
    -o /data/project/output \
    --skip-combine-gvcfs \
    --vds-path /data/project/existing.vds
```

### 10.4 With QC Filtering

```bash
# Apply QC filters before export
hvantk hgc pipeline \
    -i /data/project/gvcfs \
    -o /data/project/output \
    --apply-qc-filters \
    --min-sample-call-rate 0.90 \
    --min-variant-call-rate 0.95 \
    --generate-qc-report
```

### 10.5 Custom Temporary Directory

```bash
# Use scratch space for intermediates
hvantk hgc pipeline \
    -i /data/project/gvcfs \
    -o /data/project/output \
    --tmp-dir /scratch/tmp \
    --n-partitions 500
```

---

## 11. Appendix: Complete CLI Option Reference

```
hvantk hgc pipeline [OPTIONS]

Required Options:
  -i, --input-dir PATH           Path to directory containing input gVCF files
  -o, --output-dir PATH          Path to output directory

Stage Control:
  --skip-combine-gvcfs           Skip combining gVCF files
  --skip-vds-to-mt               Skip VDS to MatrixTable conversion
  --skip-compute-sample-qc       Skip computing sample QC metrics
  --skip-compute-variant-qc      Skip computing variant QC metrics
  --skip-export-pvcf             Skip exporting the cohort VCF

Path Overrides:
  --vds-path PATH                Path to existing VDS
  --mt-path PATH                 Path to existing MatrixTable

Configuration:
  --tmp-dir PATH                 Temporary directory for intermediate files
  --reference-genome [GRCh37|GRCh38]
                                 Reference genome build [default: GRCh38]
  --n-partitions INTEGER         Number of partitions for parallel processing
  --overwrite                    Overwrite existing output files
  --dry-run                      Show plan without executing

QC Options:
  --min-sample-call-rate FLOAT   Minimum sample call rate [default: 0.85]
  --min-variant-call-rate FLOAT  Minimum variant call rate [default: 0.85]
  --apply-qc-filters             Apply QC filters before exporting pVCF

Output Options:
  --keep-intermediates           Keep intermediate files [default: True]
  --generate-qc-report           Generate HTML QC report
  --output-prefix TEXT           Prefix for output files [default: cohort]

General:
  --help                         Show this message and exit
```

