# hvantk Architecture

## Overview

`hvantk` (Hail-based Variant Annotation Toolkit) is a modular toolkit for multi-omics variant annotation and analysis built on Apache Hail. The architecture emphasizes:

1. **Domain organization** - Separate concerns for variants, genes, proteins, and expression data
2. **Extensibility** - Protocol-based contracts for builders, streamers, and downloaders
3. **Usability** - CLI-first design with clear command structure
4. **Scalability** - Built on Hail for distributed processing of large datasets

## Project Structure

```
hvantk/
├── core/                # Shared infrastructure
│   ├── config.py       # Configuration management
│   ├── constants.py    # Shared constants
│   ├── hail_context.py # Hail session management
│   └── protocols.py    # Protocol definitions (Builder, Streamer, Downloader)
│
├── data/                # Data management utilities
│   ├── dataset.py      # Dataset handling
│   ├── file_utils.py   # File I/O utilities
│   └── data_streamer.py # Data streaming
│
├── builders/            # Data builders by domain
│   ├── variants/       # Variant-level annotations
│   │   ├── clinvar.py     # ClinVar variant annotations
│   │   ├── dbnsfp.py      # dbNSFP prediction scores
│   │   ├── gnomad.py      # gnomAD variant frequencies
│   │   └── ccr.py         # CCR constraint scores
│   │
│   ├── genes/          # Gene-level annotations
│   │   ├── ensembl.py     # Ensembl gene annotations
│   │   ├── gevir.py       # GeVIR scores
│   │   └── gnomad_metrics.py # gnomAD gene constraints
│   │
│   ├── proteins/       # Protein-level annotations
│   │   └── insider.py     # INSIDER protein-protein interactions
│   │
│   └── expression/     # Expression matrices
│       ├── ucsc.py        # UCSC Cell Browser matrices
│       └── bulk.py        # Bulk RNA-seq matrices
│
├── commands/            # CLI command implementations
│   ├── make_table_cli.py      # Table builder commands
│   ├── make_matrix_cli.py     # Matrix builder commands
│   ├── make_table_batch_cli.py  # Batch table building
│   ├── make_matrix_batch_cli.py # Batch matrix building
│   ├── catalog_cli.py         # Data catalog commands
│   ├── ucsc_downloader.py     # UCSC downloader
│   └── hgc_cli.py            # HGC joint genotyping
│
├── hgc/                 # Hail Genotype Combiner (joint genotyping)
│   ├── combine.py      # GVCF combination
│   ├── convert.py      # Format conversion (VDS ↔ MT ↔ VCF)
│   ├── qc.py           # Quality control
│   └── plot.py         # QC visualization
│
├── utils/               # Utility functions
│   ├── table_utils.py  # Table manipulation helpers
│   ├── genome.py       # Genome/contig utilities
│   └── ...
│
├── resources/           # Data catalog and schemas
│   ├── catalog.yaml    # Dataset registry
│   ├── registry/       # Dataset metadata
│   └── schemas/        # Schema definitions
│
└── tests/               # Test suite
    ├── unit/           # Unit tests
    ├── integration/    # Integration tests
    └── testdata/       # Test data fixtures
```

## Design Principles

### 1. Domain Separation

Data is organized by biological domain rather than technical implementation:

- **Variants** (`builders/variants/`) - Variant-level annotations keyed by `(locus, alleles)`
- **Genes** (`builders/genes/`) - Gene-level annotations keyed by `gene_id`
- **Proteins** (`builders/proteins/`) - Protein-level annotations keyed by `protein_id` or `interval`
- **Expression** (`builders/expression/`) - Expression matrices with rows=genes, columns=samples/cells

This organization makes it intuitive to find and add new data sources.

### 2. Protocol-Based Extensibility

Three core protocols define how components interact:

#### Builder Protocol
Converts raw data files → Hail Tables/MatrixTables

```python
from hvantk.core.protocols import Builder

class MyBuilder:
    def build(self, input_path: str, **params) -> hl.Table:
        """Convert raw file to Hail Table"""
        
    def validate_schema(self, ht: hl.Table) -> bool:
        """Validate output schema"""
        
    def get_metadata(self) -> Dict[str, Any]:
        """Return builder metadata"""
```

#### Streamer Protocol
Transforms Hail data structures (filter, join, aggregate)

```python
from hvantk.core.protocols import Streamer

class MyStreamer:
    def transform(self, input_data: hl.Table, **params) -> hl.Table:
        """Transform input data"""
        
    def validate_input(self, input_data: hl.Table) -> bool:
        """Validate input schema"""
        
    def get_metadata(self) -> Dict[str, Any]:
        """Return streamer metadata"""
```

#### Downloader Protocol
Fetches external datasets with verification

```python
from hvantk.core.protocols import Downloader
from pathlib import Path

class MyDownloader:
    def download(self, dataset_id: str, output_dir: Path, **params) -> Path:
        """Download dataset"""
        
    def verify_checksum(self, file_path: Path, expected: str) -> bool:
        """Verify file integrity"""
        
    def get_metadata(self, dataset_id: str) -> Dict[str, Any]:
        """Return dataset metadata"""
```

### 3. CLI-First Design

The primary interface is a well-structured CLI with domain-specific commands:

```bash
# Download data
hvantk ucsc-downloader --dataset adultPancreas --output-dir data/

# Build individual tables
hvantk mktable clinvar --raw-input clinvar.vcf.bgz --output-ht clinvar.ht
hvantk mktable ensembl-gene --raw-input biomart.tsv --output-ht ensembl.ht

# Build matrices
hvantk mkmatrix ucsc -e expr.tsv.bgz -m meta.tsv -o ucsc.mt

# Batch processing via recipes
hvantk mktable-batch --recipe tables.json
hvantk mkmatrix-batch --recipe matrices.json

# Joint genotyping (HGC)
hvantk hgc gvcf-combine -g /data/gvcfs -o cohort.vds
hvantk hgc compute-qc -i cohort.mt -o cohort_qc.mt
```

### 4. Data Flow Patterns

#### Pattern 1: Builder → Table
```
Raw File (VCF/TSV/BED) → Builder → Hail Table → Disk (.ht)
```

Example:
```python
from hvantk.builders.variants.clinvar import create_clinvar_tb

ht = create_clinvar_tb(
    input_path="clinvar.vcf.bgz",
    output_path="clinvar.ht",
    reference_genome="GRCh38"
)
```

#### Pattern 2: Batch Building via Recipes
```
Recipe JSON → Parser → [Builder 1, Builder 2, ...] → Multiple Tables
```

Example recipe:
```json
{
  "tables": [
    {
      "name": "clinvar",
      "input": "/data/clinvar.vcf.bgz",
      "output": "/out/clinvar.ht",
      "params": {"reference_genome": "GRCh38"}
    }
  ]
}
```

#### Pattern 3: Multi-Omics Integration
```
Variant Table + Gene Table + Expression Matrix → Integrated Analysis
```

Example:
```python
# Load different data types
variants = hl.read_table('clinvar.ht')
genes = hl.read_table('ensembl.ht')
expression = hl.read_matrix_table('ucsc.mt')

# Join for integrated analysis
annotated = variants.annotate(
    gene=genes[variants.gene_id],
    expression=expression[variants.gene_id, :].expression.collect()
)
```

## Module Details

### Core Module (`core/`)

**Purpose**: Shared infrastructure used by all other modules

**Key Components**:
- `config.py` - Configuration management, context settings
- `constants.py` - Shared constants (e.g., Ensembl field definitions)
- `hail_context.py` - Hail session initialization and management
- `protocols.py` - Protocol definitions for extensibility

**Design principle**: No domain logic, only infrastructure

### Data Module (`data/`)

**Purpose**: Data management utilities

**Key Components**:
- `dataset.py` - Dataset handling and metadata
- `file_utils.py` - File I/O utilities (download, checksum, compression)
- `data_streamer.py` - Data streaming and transformation helpers

### Builders Module (`builders/`)

**Purpose**: Convert raw data files into Hail Tables/MatrixTables

**Organization by domain**:

#### `builders/variants/`
- **ClinVar** - Variant clinical significance (VCF → Table)
- **dbNSFP** - Missense variant prediction scores (TSV → Table)
- **gnomAD** - Population variant frequencies (VCF → Table)
- **CCR** - Coding-constrained regions (BED → Table)

**Schema**: Tables keyed by `(locus, alleles)`

#### `builders/genes/`
- **Ensembl** - Gene annotations from Biomart (TSV → Table)
- **GeVIR** - Gene-level viability scores (TSV → Table)
- **gnomAD Metrics** - Gene constraint metrics (TSV → Table)

**Schema**: Tables keyed by `gene_id`

#### `builders/proteins/`
- **INSIDER** - Protein-protein interaction sites (BED → Table)

**Schema**: Tables keyed by `interval` or `protein_id`

#### `builders/expression/`
- **UCSC** - Single-cell RNA-seq from UCSC Cell Browser (TSV → MatrixTable)
- **Bulk** - Bulk RNA-seq expression (TSV → MatrixTable)

**Schema**: MatrixTables with rows=genes, columns=samples/cells

### Commands Module (`commands/`)

**Purpose**: CLI command implementations

**Key files**:
- `make_table_cli.py` - Commands for building individual tables
- `make_matrix_cli.py` - Commands for building matrices
- `make_table_batch_cli.py` - Batch table building from recipes
- `make_matrix_batch_cli.py` - Batch matrix building from recipes
- `catalog_cli.py` - Data catalog operations
- `hgc_cli.py` - HGC joint genotyping commands

### HGC Module (`hgc/`)

**Purpose**: High-performance joint genotyping workflows

**Features**:
- GVCF combination at scale
- Format conversion (VDS ↔ MatrixTable ↔ VCF)
- Comprehensive QC metrics and visualization
- Optimized for large cohorts (1000s of samples)

**Status**: Feature-complete, well-established module

### Resources Module (`resources/`)

**Purpose**: Data catalog and schema definitions

**Contents**:
- `catalog.yaml` - Central dataset registry
- `registry/` - Per-domain dataset metadata (genomics, transcriptomics, etc.)
- `schemas/` - Schema definitions for validation

## Testing Strategy

Tests are organized to mirror the module structure:

```
tests/
├── unit/                # Unit tests for individual components
│   ├── builders/       # Builder tests
│   ├── test_core.py   # Core utilities
│   └── test_commands.py # CLI commands
│
├── integration/         # Integration tests for workflows
│   └── test_workflows.py
│
└── testdata/            # Test fixtures
    ├── raw/            # Sample raw data files
    └── expected/       # Expected output schemas
```

## Extension Points

### Adding a New Data Source

1. **Choose the appropriate domain** (`variants/`, `genes/`, `proteins/`, `expression/`)

2. **Implement a builder** following the `Builder` protocol:
   ```python
   # hvantk/builders/variants/my_source.py
   from hvantk.core.protocols import Builder
   import hail as hl
   
   def create_my_source_tb(input_path: str, output_path: str, **kwargs) -> hl.Table:
       # Implementation
       pass
   ```

3. **Add a CLI command**:
   ```python
   # hvantk/commands/make_table_cli.py
   @mktable_group.command("my-source")
   def mktable_my_source(...):
       create_my_source_tb(...)
   ```

4. **Add tests**:
   ```python
   # hvantk/tests/unit/builders/test_variants.py
   def test_create_my_source_tb():
       # Test implementation
       pass
   ```

5. **Update documentation** in README.md and USAGE.md

### Adding a New Transformation

1. **Implement a streamer** following the `Streamer` protocol
2. **Add to the pipeline** (for batch processing support)
3. **Document** the transformation parameters

## Dependencies

- **Apache Hail** - Distributed data processing framework
- **gnomAD** - Utilities for gnomAD data
- **Click** - CLI framework
- **Pandas** - Data manipulation
- **PyYAML** - YAML recipe support
- **Matplotlib/Seaborn/Plotly** - Visualization (optional)

## Performance Considerations

- **Partitioning**: Builders use appropriate partitioning for Hail operations
- **Checkpointing**: Large intermediate results are checkpointed
- **Memory**: HGC module optimized for memory-efficient large cohort processing
- **Caching**: Hail's lazy evaluation allows for optimization

## Future Directions

1. **More data sources**: Add support for additional variant/gene databases
2. **Streamers library**: Build reusable transformation components
3. **Pipeline DAGs**: Support for complex multi-step workflows
4. **Cloud integration**: Better support for S3/GCS data sources
5. **Web API**: Optional REST API for programmatic access

## References

- [Hail Documentation](https://hail.is/docs/0.2/)
- [gnomAD Browser](https://gnomad.broadinstitute.org/)
- [UCSC Cell Browser](https://cells.ucsc.edu/)
