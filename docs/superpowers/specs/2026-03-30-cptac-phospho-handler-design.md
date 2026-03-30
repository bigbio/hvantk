# CPTAC Phospho Handler — Design Spec

**Date:** 2026-03-30
**Status:** Approved
**Depends on:** PeptideAtlas phospho handler (PTM_OUTPUT_COLUMNS 13-col format)

## Summary

Add CPTAC (Clinical Proteomic Tumor Analysis Consortium) as a phospho site data source in the PTM pipeline, using the `cptac` Python package to access phosphoproteomics data from 7+ cancer types. Produces two outputs: (1) a site-level summary TSV for PTM pipeline integration with genomic coordinate mapping, and (2) a phospho MatrixTable (sites × samples) with per-sample log2 intensities and clinical metadata.

## Data Source

The `cptac` Python package (PyPI: `pip install cptac`) provides programmatic access to CPTAC datasets. Each cancer type exposes `get_phosphoproteomics()` returning a pandas DataFrame with:
- Rows: samples (~100 tumor + ~30 normal per cancer type)
- Columns: multi-index (Gene, Site, Peptide, Database_ID)
- Values: log2 intensity ratios (TMT/iTRAQ)

Clinical metadata via `get_clinical()` provides sample_type (tumor/normal), stage, subtype, age, etc.

### Available Cancer Types

brca, ccrcc, colon, endometrial, gbm, hnscc, lscc, luad, ov, pdac, ucec

### Site Nomenclature

```
Gene    Database_ID    Phosphosite    Peptide
EIF4EBP1  NP_004086.1  S65           RVpSGGEELGS
TP53      NP_000537.3  S315          SPQPKKKPLDGEpS
MAPK1     NP_002736.3  T185_Y187    VADPDHDHTGFLpTEpYVATR
```

**Parsing rules:**
- Single site: `S65` → amino_acid="S", position=65
- Multi-site: `T185_Y187` → two separate entries, one per site
- Only S, T, Y residues accepted (consistent with PeptideAtlas)

## Output Schema

### Output 1: PTM Pipeline Intermediate TSV

Same 13-column mapped format as UniProt/PeptideAtlas, extended with `mean_intensity`:

| Column | Value for CPTAC |
|--------|----------------|
| chrom | (after GTF mapping) |
| codon_start | (after GTF mapping) |
| codon_end | (after GTF mapping) |
| strand | (after GTF mapping) |
| uniprot_id | (empty — resolved via gene_mane) |
| gene_symbol | Gene name from cptac |
| residue_pos | Parsed from site column |
| amino_acid | S, T, or Y |
| ptm_type | Phosphoserine/Phosphothreonine/Phosphotyrosine |
| ptm_category | phosphorylation |
| source_db | CPTAC |
| evidence_type | mass_spectrometry |
| n_observations | Number of samples with non-NaN intensity |

The intermediate TSV (pre-mapping) includes the standard 10 columns plus `n_observations`, `source_db`, `evidence_type`, and adds `mean_intensity` and `cancer_type`.

### Output 2: Phospho MatrixTable

- **Rows:** phospho sites keyed by (gene_symbol, residue_pos, amino_acid)
- **Row annotations:** ptm_type, ptm_category, database_id, peptide
- **Columns:** samples keyed by sample_id
- **Column annotations:** cancer_type, sample_type (Tumor/Normal), clinical metadata from `get_clinical()`
- **Entry:** log2 intensity (Float64)

## Data Flow

```
cptac Python package
  → CPTACPhosphoDataset(cancer_type="brca")
  → dataset.download(output_dir)
    1. cptac.Brca().get_phosphoproteomics() → DataFrame
    2. cptac.Brca().get_clinical() → DataFrame
    3. Parse site columns → (gene, aa, pos) tuples
    4. Aggregate per site: n_samples_detected, mean_intensity
    5. Write intermediate TSV (compatible with PTM pipeline mapper)
    6. Write phospho matrix CSV (sites × samples for MatrixTable)
    7. Write clinical metadata CSV

  PTM pipeline integration:
    hvantk ptm build --cptac-tsv intermediate.tsv
    → GTF mapper (gene_mane fallback) → 13-col mapped TSV → Hail Table

  MatrixTable path:
    hvantk mkmatrix cptac-phospho -e matrix.csv -m metadata.csv -o phospho.mt
```

## Transcript Resolution

CPTAC provides gene symbols but no Ensembl cross-references. The intermediate TSV sets `ensembl_xrefs=""`, and the existing GTF mapper resolves via the `gene_mane` fallback (gene symbol → MANE Select transcript), the same strategy used for PeptideAtlas.

## CLI Commands

### Download

```bash
# Single cancer type
hvantk download cptac-phospho --cancer-type brca -o data/cptac/

# All cancer types (pan-cancer)
hvantk download cptac-phospho --all -o data/cptac/

# List available cancer types
hvantk download cptac-phospho --list-cancers
```

Options:
- `--cancer-type` — One of: brca, ccrcc, colon, endometrial, gbm, hnscc, lscc, luad, ov, pdac, ucec
- `--all` — Download and merge all cancer types
- `--list-cancers` — Print available cancer types and exit
- `--output-dir` (required)
- `--overwrite`

Output files per cancer type:
- `cptac-phospho-{cancer_type}.tsv` — intermediate TSV for PTM pipeline
- `cptac-phospho-{cancer_type}-matrix.csv` — sites × samples matrix
- `cptac-phospho-{cancer_type}-metadata.csv` — sample clinical metadata

With `--all`: additional `cptac-phospho-pancancer.tsv` combining all cancer types with `cancer_type` column.

### PTM Build Integration

```bash
hvantk ptm build -o data/ptm/ --output-ht ptm.ht \
    --cptac-tsv data/cptac/cptac-phospho-brca.tsv
```

### MatrixTable Build

```bash
hvantk mkmatrix cptac-phospho \
    -e data/cptac/cptac-phospho-brca-matrix.csv \
    -m data/cptac/cptac-phospho-brca-metadata.csv \
    -o data/cptac/brca_phospho.mt
```

## Files to Create

| File | Responsibility |
|------|----------------|
| `hvantk/datasets/cptac_phospho_datasets.py` | Dataset class: wrap cptac package, parse sites, aggregate, write TSV + matrix |
| `hvantk/commands/cptac_phospho_downloader.py` | CLI `hvantk download cptac-phospho` |
| `hvantk/tests/test_cptac_phospho.py` | Unit tests with mock DataFrames |

## Files to Modify

| File | Change |
|------|--------|
| `hvantk/ptm/constants.py` | Add CPTAC constants (cancer types, site parsing patterns) |
| `hvantk/ptm/pipeline.py` | Add `cptac_tsv` to PTMBuildConfig, concatenation in ptm_build_pipeline() |
| `hvantk/commands/ptm_cli.py` | Add `--cptac-tsv` option to ptm build |
| `hvantk/commands/download_cli.py` | Register cptac-phospho downloader |
| `hvantk/tables/matrix_builders.py` | Add `build_cptac_phospho_mt()` |
| `hvantk/tables/registry.py` | Register "cptac-phospho" in MATRIX_BUILDERS |
| `hvantk/commands/make_matrix_cli.py` | Add `hvantk mkmatrix cptac-phospho` command |

## Clinical Metadata Columns

Standard clinical columns captured from `get_clinical()`:
- `sample_type` — Tumor / Normal
- `age` — Patient age
- `gender` — Patient gender
- `tumor_stage` — TNM staging (when available)
- Additional cancer-type-specific columns preserved as-is

## Dependency

The `cptac` Python package is an optional dependency (lazy-imported). If not installed, the downloader raises a clear error with install instructions.

## Backward Compatibility

- CPTAC rows in the PTM pipeline: `source_db="CPTAC"`, `evidence_type="mass_spectrometry"`, `n_observations=<n_samples>`
- Existing UniProt and PeptideAtlas rows unaffected
- The existing `mkmatrix cptac` (gene-level expression) is untouched; the new `mkmatrix cptac-phospho` is a separate command
