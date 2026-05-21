# hvantk Architecture

## Overview

`hvantk` (Hail-based Variant Annotation Toolkit) is a modular toolkit for multi-omics variant annotation and analysis built on Hail. The architecture emphasizes:

1. **Domain organization** - Separate concerns for variants, genes, proteins, and expression data
2. **Extensibility** - Protocol-based contracts for builders, streamers, and downloaders
3. **Usability** - CLI-first design with clear command structure
4. **Scalability** - Built on Hail for distributed processing of large datasets

## Architecture Diagram

![hvantk workflow architecture](images/hvantk-architecture.svg)

**Figure 1.** *hvantk workflow architecture for scalable multi-omics variant annotation and analysis. External variant, gene, and expression databases are acquired through built-in downloaders and converted to domain-organized Hail Tables and MatrixTables via the builder framework. User cohort data (GVCFs, MatrixTables, gene sets, phenotypes) feeds directly into analysis pipelines. Four specialized pipelines — HGC (joint genotyping and quality control), Ancestry (PCA-based population inference), EnrichEx (gene-set burden and overlap testing), and PS-ROC (pathogenicity score evaluation) — produce annotated tables, HTML reports with embedded plots, and statistical results. All operations are distributed via Hail on Apache Spark, accessible through the `hvantk` CLI and Python API.*

## Project Structure

```text
hvantk/
├── __main__.py            # Main CLI entry point
│
├── core/                  # L1: Core infrastructure
│   ├── config.py          # Configuration management
│   ├── constants.py       # Shared constants
│   ├── protocols.py       # Protocol definitions (Builder, Streamer, Downloader)
│   ├── builders/          # Generic builder helpers
│   │   └── table.py       # _create_table_base, _cleanup_temp_file, etc.
│   ├── io/                # Artifact loader (load/save Hail Tables, AnnData, etc.)
│   ├── models/            # Domain model types
│   │   ├── annotation_table.py  # AnnotationTable artifact
│   │   ├── expression_matrix.py # ExpressionMatrix artifact
│   │   ├── gene_set.py          # GeneSet artifact
│   │   ├── artifact.py          # Artifact base + type registry
│   │   ├── backends.py          # AlgorithmMeta, Backend, @algorithm decorator
│   │   ├── build_context.py     # BuildContext passed to plugin builders
│   │   ├── anndata_utils.py     # AnnData helpers (build_anndata_metadata, etc.)
│   │   ├── metadata.py          # Metadata structs and source descriptions
│   │   └── provenance.py        # Source-fingerprint provenance stamping
│   ├── plugin/            # Plugin system
│   │   ├── api.py         # PluginSpec, DatasetSpec, DriftProbeError
│   │   ├── loader.py      # Plugin discovery (filesystem + entry points)
│   │   ├── registry.py    # TABLE_BUILDERS / MATRIX_BUILDERS registry
│   │   ├── run_builder.py # run_builder_for_spec() — validates artifact type
│   │   └── drift_runner.py# Drift probe execution
│   ├── streamers/         # Data streamers (base, gene_disease, etc.)
│   └── utils/             # Cross-cutting utilities
│       ├── hail_context.py  # Idempotent Hail init
│       ├── bgzf.py          # BGZF utilities
│       ├── file_utils.py    # File I/O helpers
│       ├── gene_sets.py     # Gene set utilities
│       ├── genome.py        # Genome/contig utilities
│       ├── table_utils.py   # Hail Table manipulation helpers
│       ├── writers.py       # HailTableWriter
│       └── ...              # Other shared utilities
│
├── algorithms/            # L4-L5: Analysis pipelines
│   ├── annotation/        # Variant annotation pipeline
│   ├── ancestry/          # Population ancestry inference
│   ├── enrichex/          # Gene set enrichment analysis
│   ├── expression/        # Expression data processing
│   ├── hgc/               # Joint genotyping (HGC) pipeline
│   ├── psroc/             # Pathogenicity score evaluation
│   ├── ptm/               # Post-translational modification analysis
│   ├── qtlcascade/        # QTL cascade analysis
│   ├── statistics/        # Statistical utilities
│   ├── training_sets/     # Training set construction
│   └── visualization/     # Shared visualization helpers
│
├── skills/                # L2-L3: Per-provider data plugins
│   ├── _conventions/      # Shared plugin contract (SKILL.md)
│   ├── _hooks/            # Plugin lifecycle hooks
│   ├── clingen/           # ClinGen gene-disease validity
│   ├── clinvar/           # ClinVar variant annotations
│   ├── cptac/             # CPTAC proteomics (expression/, phospho/)
│   ├── expression_atlas/  # Expression Atlas bulk RNA-seq
│   ├── gencc/             # GenCC gene-disease assertions
│   ├── gtex_eqtl/         # GTEx eQTL data
│   ├── gwas_catalog/      # GWAS Catalog
│   ├── hgnc/              # HGNC gene nomenclature
│   ├── insider/           # INSIDER protein-protein interaction sites
│   ├── msigdb/            # MSigDB gene sets
│   ├── peptideatlas/      # PeptideAtlas proteomics
│   ├── ucsc_cellbrowser/  # UCSC Cell Browser single-cell RNA-seq
│   └── uniprot_ptm/       # UniProt PTM annotations
│
├── tools/                 # CLI command implementations (replaces legacy commands/)
│   ├── ancestry/          # Ancestry CLI subcommands
│   ├── annotation/        # Annotation CLI subcommands
│   ├── build/             # mktable / mkmatrix / batch build commands
│   ├── enrichex/          # EnrichEx CLI subcommands
│   ├── expression/        # Expression analysis commands
│   ├── genesets/          # Gene set extraction/preparation
│   ├── hgc/               # HGC CLI subcommands
│   ├── infra/             # Installation check, BGZF validation, utils
│   ├── plugins/           # hvantk plugins / hvantk drift commands
│   ├── ptm/               # PTM CLI subcommands
│   └── qtl/               # QTL CLI subcommands
│
├── resources/             # Data catalog and schemas
│   ├── registry/          # Surviving legacy per-domain dataset metadata (genomics only)
│   ├── schemas/           # JSON schema definitions
│   └── unified_registry.py# Aggregates per-plugin catalog/datasets.json + legacy registry
│
└── tests/                 # Test suite
    ├── conftest.py        # Pytest fixtures (hail_session, etc.)
    ├── testdata/          # Test data fixtures
    ├── hgc/               # HGC tests
    ├── ancestry/          # Ancestry tests
    ├── psroc/             # PSROC tests
    ├── enrichex/          # EnrichEx tests
    └── test_*.py          # Unit and integration tests
```

## Design Principles

### 1. Domain Separation

The codebase is organized by function and biological domain:

**Data Builders** (`skills/<provider>/builder.py`):
- Each plugin under `hvantk/skills/` owns its builder. Builders return `AnnotationTable`, `ExpressionMatrix`, or `GeneSet` artifacts.
- Generic helpers live in `hvantk/core/builders/table.py` (`_create_table_base`, etc.).

**Analysis Pipelines** (separate modules):
- `hgc/` - Joint genotyping and cohort analysis
- `ancestry/` - Population ancestry inference
- `psroc/` - Pathogenicity score evaluation
- `enrichex/` - Gene set enrichment analysis
- `ptm/` - Post-translational modification variant classification. Includes
  `constraint.py` + helpers for the stratified AF-depletion analysis
  (`hvantk ptm constraint`), a tissue/cell-type-aware complement to
  `landscape` and `population`.

**Data Product Keying**:
- **Variants** - Keyed by `(locus, alleles)`
- **Genes** - Keyed by `gene_id`
- **Proteins** - Keyed by `protein_id` or `interval`
- **Expression** - MatrixTables with rows=genes, columns=samples/cells

### 2. Protocol-Based Extensibility

Three core protocols define how components interact:

#### Builder Protocol (legacy — pre-Phase B)

> **Deprecated.** This section describes the pre-Phase B builder shape
> (`(input_path, output_path) -> hl.Table`). New plugin authors should
> use the Phase B contract documented in the "Adding a New Data Source"
> section below: `(parsed_input, ctx: BuildContext, **params) -> Artifact`.
> Legacy functions in `hvantk/core/builders/table.py` are retained for
> backward compatibility with existing recipes; new builders live in
> `hvantk/skills/<plugin>/builder.py` and return artifact instances.

Converts raw data files → Hail Tables/MatrixTables

Builders follow a functional pattern using `_create_table_base()` to eliminate boilerplate:

```python
import hail as hl
from hvantk.core.builders.table import _create_table_base

def create_my_source_tb(input_path: str, output_path: str, **kwargs) -> hl.Table:
    """Build a Hail Table from MySource data.

    Assumes imported records contain `locus` and `alleles` fields.
    """
    return _create_table_base(
        source_name="MySource",
        input_path=input_path,
        output_path=output_path,
        import_func=lambda: hl.import_table(input_path, ...),
        transform_func=lambda ht: ht.key_by(ht.locus, ht.alleles),
        overwrite=kwargs.get('overwrite', False),
        export_tsv=kwargs.get('export_tsv', False),
    )
```

#### Streamer Protocol
Transforms Hail data structures (filter, join, aggregate)

Streamers extend `HailDataStreamer` from `hvantk/core/streamers/base.py`:

```python
from typing import Iterator
import hail as hl
from hvantk.core.streamers.base import HailDataStreamer

class MySourceStreamer(HailDataStreamer):
    def __init__(self, table_path: str, chunk_size: int = 10000):
        super().__init__("MySourceStreamer", chunk_size=chunk_size)
        self.table_path = table_path

    def setup(self) -> None:
        super().setup()
        self._table = hl.read_table(self.table_path)

    def stream(self) -> Iterator[hl.Table]:
        # Yield chunks of data
        ...
```

#### Downloader Protocol
Fetches external datasets with verification

Downloaders use dataset dataclasses with a `download()` method:

```python
from dataclasses import dataclass
from pathlib import Path

@dataclass
class MyDataset:
    url: str
    output_dir: Path

    def download(self, overwrite: bool = False) -> Path:
        """Download and verify dataset."""
        ...

    @classmethod
    def latest(cls, output_dir: Path) -> "MyDataset":
        """Create instance for the latest available version."""
        ...
```

### 3. CLI-First Design

The primary interface is a well-structured CLI with domain-specific commands:

```bash
# Download data
hvantk download ucsc --dataset adultPancreas --output-dir data/

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
# Via the plugin system (recommended)
from hvantk.core.plugin.run_builder import run_builder_for_spec

artifact = run_builder_for_spec("clinvar:variants", input_path="clinvar.vcf.bgz", output_path="clinvar.ht")
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
- `utils/hail_context.py` - Hail session initialization and management
- `protocols.py` - Protocol definitions for extensibility
- `models/` - Domain artifact types (`AnnotationTable`, `ExpressionMatrix`, `GeneSet`)
- `plugin/` - Plugin schema (`api.py`), discovery (`loader.py`), and builder dispatch (`run_builder.py`)

**Design principle**: No domain logic, only infrastructure

### Skills Module (`skills/`)

**Purpose**: Per-provider data plugins. Each provider folder contains `plugin.yaml`, `builder.py`, `cli.py`, `drift_probe.py`, `SKILL.md`, `catalog/datasets.json`, and `tests/`. Multi-dataset providers (e.g., `cptac/`) have one sub-folder per dataset.

**Current providers**: `clingen`, `clinvar`, `cptac`, `expression_atlas`, `gencc`, `gtex_eqtl`, `gwas_catalog`, `hgnc`, `insider`, `msigdb`, `peptideatlas`, `ucsc_cellbrowser`, `uniprot_ptm`.

**Builder outputs**:
- Variant / gene tables keyed by `(locus, alleles)` or `gene_id` → `AnnotationTable`
- Expression matrices rows=genes, columns=samples/cells → `ExpressionMatrix`
- Gene set collections → `GeneSet`

### Tools Module (`tools/`)

**Purpose**: Top-level CLI command implementations (replaces the legacy `commands/` directory)

**Key sub-packages**:
- `build/` - `mktable`, `mkmatrix`, `mktable-batch`, `mkmatrix-batch` commands
- `plugins/` - `hvantk plugins list/show/reload` and `hvantk drift` commands
- `hgc/` - HGC joint genotyping subcommands (combine, convert, QC, pipeline)
- `ancestry/`, `enrichex/`, `ptm/`, `qtl/` - Per-pipeline CLI subcommands
- `infra/` - Installation check, BGZF validation, utils

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
- `registry/` - Surviving legacy per-domain dataset metadata (genomics only; transcriptomics / proteomics / epigenomics moved into per-plugin `hvantk/skills/<provider>/catalog/datasets.json`)
- `unified_registry.py` - `HvantkRegistry` aggregator surfaced via `hvantk catalog {list,show,stats,search}`
- `schemas/` - Schema definitions for validation

## Testing Strategy

Tests are organized to mirror the module structure:

```
hvantk/tests/
├── conftest.py        # Pytest fixtures (hail_session, etc.)
├── test_*.py          # Unit and integration tests
├── hgc/               # HGC module tests
├── ancestry/          # Ancestry module tests
├── psroc/             # PSROC module tests
├── enrichex/          # EnrichEx module tests
└── testdata/          # Test data fixtures
```

## Extension Points

### Adding a New Data Source

1. **Create a new plugin folder** under `hvantk/skills/<provider>/` with `plugin.yaml`, `builder.py`, `cli.py`, `drift_probe.py`, `SKILL.md`, `catalog/datasets.json`, and `tests/`. See `hvantk/skills/_conventions/SKILL.md` for the full contract.

2. **Implement the builder** in `hvantk/skills/<provider>/builder.py`:
   ```python
   import hail as hl
   from hvantk.core.builders.table import _create_table_base
   from hvantk.core.models.build_context import BuildContext

   def create_my_source_tb(parsed_input, ctx: BuildContext, **kwargs):
       """
       Create an AnnotationTable from my data source.
       Returns an AnnotationTable artifact.
       """
       return _create_table_base(
           source_name="MySource",
           input_path=parsed_input.path,
           output_path=ctx.output_path,
           import_func=lambda: hl.import_table(parsed_input.path, ...),
           transform_func=lambda ht: ht.key_by(ht.locus, ht.alleles),
           overwrite=kwargs.get('overwrite', False),
       )
   ```

3. **Add a CLI command** in `hvantk/skills/<provider>/cli.py` and declare it in `plugin.yaml`:
   ```yaml
   cli:
     - command: my-source-download
       module: hvantk.skills.my_source.cli
       function: download_cmd
   ```

4. **Add tests** in `hvantk/skills/<provider>/tests/`:
   ```python
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

- **Hail** - Distributed data processing framework
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

## References

- [Hail Documentation](https://hail.is/docs/0.2/)
- [gnomAD Browser](https://gnomad.broadinstitute.org/)
- [UCSC Cell Browser](https://cells.ucsc.edu/)
