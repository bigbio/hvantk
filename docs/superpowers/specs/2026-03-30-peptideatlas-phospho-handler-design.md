# PeptideAtlas Phospho Handler — Design Spec

**Date:** 2026-03-30
**Status:** Draft
**Source:** https://peptideatlas.org/builds/human/phospho/

## Summary

Add PeptideAtlas human phospho build as a data source in the existing PTM pipeline. Downloads the TSV dump zip, parses phospho site positions with observation counts, maps protein positions to genomic coordinates using the existing GTF-based mapper, and outputs a 13-column mapped TSV compatible with the current pipeline.

## Output Schema (13 columns)

Extends the current 12-column `PTM_OUTPUT_COLUMNS` with `n_observations`:

| Column          | Type   | Description                                      |
|-----------------|--------|--------------------------------------------------|
| chrom           | str    | Chromosome (e.g., "chr17")                       |
| codon_start     | int    | Lowest genomic position of 3-bp codon            |
| codon_end       | int    | Highest genomic position of 3-bp codon           |
| strand          | str    | "+" or "-"                                       |
| uniprot_id      | str    | UniProt accession                                |
| gene_symbol     | str    | Gene symbol                                      |
| residue_pos     | int    | 1-based protein residue position                 |
| amino_acid      | str    | Amino acid at PTM site (S, T, Y for phospho)     |
| ptm_type        | str    | PTM description (e.g., "Phosphoserine")          |
| ptm_category    | str    | Normalized category ("phosphorylation")          |
| source_db       | str    | "UniProt" or "PeptideAtlas"                      |
| evidence_type   | str    | "curated" (UniProt) or "mass_spectrometry" (PA)  |
| n_observations  | int    | PSM/observation count (0 for UniProt rows)       |

## Data Flow

```
PeptideAtlas TSV zip (523MB)
  → download & extract
  → parse tables: biosequence, peptide_mapping, modified_peptide_instance, peptide_instance
  → aggregate: per (accession, residue_position) → sum(n_observations)
  → intermediate TSV: accession, gene_symbol, position, description, amino_acid,
    ensembl_xrefs, sequence_length, n_observations
  → existing GTF mapper (resolve_transcript → map_protein_sites)
  → 13-column mapped TSV
  → Hail Table via updated create_ptm_sites_tb
```

## PeptideAtlas TSV Parsing Strategy

1. Download zip from `http://www.peptideatlas.org/builds/{BUILD_DATE}/atlas_build_{BUILD_ID}.tsv.zip`
2. Extract and discover table files by name pattern
3. Identify tables: `biosequence`, `peptide_mapping`, `modified_peptide_instance`, `peptide_instance`
4. Read headers to discover exact column names
5. Join logic:
   - `biosequence` → UniProt accession, gene name, sequence
   - `peptide_mapping` → maps peptide to protein position (start_in_biosequence)
   - `modified_peptide_instance` → phospho modification offset within peptide
   - `peptide_instance` → n_observations, n_samples
   - Site position = `start_in_biosequence + modification_offset_within_peptide`
6. Aggregate by (accession, site_position) → sum of n_observations
7. Fallback: if table structure differs from expected, log discovered tables/columns and raise a clear error

## Transcript Resolution Note

PeptideAtlas biosequence table contains UniProt accessions but may not carry Ensembl cross-references. The existing transcript resolver has a 3-strategy cascade: `xref_mane` → `xref_any` → `gene_mane`. For PeptideAtlas proteins, the resolver will primarily use the `gene_mane` fallback (gene symbol → MANE Select transcript). If Ensembl xrefs are available in the biosequence table, they'll be used; otherwise the gene name path handles it.

## Build Versioning

URL pattern: `http://www.peptideatlas.org/builds/{BUILD_DATE}/atlas_build_{BUILD_ID}.tsv.zip`

The dataset class takes `build_date` and `build_id` parameters. Defaults point to the latest build (2025-12, build 606). Constants stored in `hvantk/ptm/constants.py`.

## Files to Create

### 1. `hvantk/datasets/peptideatlas_phospho_datasets.py`

`PeptideAtlasPhosphoDataset` class following the `UniProtPTMDataset` pattern:
- `from_latest()` / `from_build(build_date, build_id)` class methods
- `download(output_dir, overwrite)` — downloads zip, extracts, parses tables, aggregates phospho sites, writes intermediate TSV
- `get_metadata()` — returns dataset metadata dict
- Internal helpers: `_parse_biosequences()`, `_parse_peptide_mappings()`, `_parse_phospho_sites()`, `_aggregate_sites()`

Output: intermediate TSV with columns matching UniProt's 7-column format + `n_observations`.

### 2. `hvantk/commands/peptideatlas_phospho_downloader.py`

CLI command: `hvantk download peptideatlas-phospho`

Options:
- `--output-dir` (required)
- `--build-date` (default: latest)
- `--build-id` (default: latest)
- `--overwrite`

### 3. `hvantk/tests/test_peptideatlas_phospho.py`

Unit tests for:
- TSV parsing with mock data (small synthetic tables)
- Phospho site aggregation logic
- Edge cases: multiple peptides covering same site, missing accessions
- Dataset class construction and metadata

## Files to Modify

### 4. `hvantk/ptm/constants.py`

Add:
```python
# PeptideAtlas Phospho Build
PEPTIDEATLAS_PHOSPHO_BASE_URL = "http://www.peptideatlas.org/builds"
PEPTIDEATLAS_LATEST_BUILD_DATE = "202512"
PEPTIDEATLAS_LATEST_BUILD_ID = "606"
```

Update:
```python
PTM_OUTPUT_COLUMNS = [
    "chrom", "codon_start", "codon_end", "strand", "uniprot_id",
    "gene_symbol", "residue_pos", "amino_acid", "ptm_type",
    "ptm_category", "source_db", "evidence_type", "n_observations",
]
```

### 5. `hvantk/ptm/pipeline.py`

- Update `map_ptm_sites()` to pass through `source_db`, `evidence_type`, and `n_observations` from input to output
- Add `download_peptideatlas_phospho()` helper
- Update `ptm_build_pipeline()` to optionally accept a PeptideAtlas TSV, map it, and concatenate with UniProt mapped output before building the Hail Table

### 6. `hvantk/commands/ptm_cli.py`

- Add `--peptideatlas-tsv` option to `ptm build` command
- When supplied, pipeline maps PeptideAtlas sites alongside UniProt sites

### 7. `hvantk/commands/download_cli.py` (or equivalent registration point)

Register `peptideatlas_phospho_downloader` in `download_group`.

### 8. `hvantk/tables/table_builders.py`

Update `create_ptm_sites_tb` to:
- Import `n_observations` as `Int32` (default 0)
- Include in Hail Table schema

## Backward Compatibility

- UniProt rows get `evidence_type="curated"`, `n_observations=0`
- Existing downstream analysis (annotate, landscape, population) is unaffected — they key on `locus` and `ptm_category`, not on the new fields
- The 13-column format is a strict superset of the 12-column format

## CLI Usage

```bash
# Download PeptideAtlas phospho build
hvantk download peptideatlas-phospho -o data/peptideatlas/

# Build PTM table from PeptideAtlas only
hvantk ptm build -o data/ptm/ --output-ht ptm.ht \
    --peptideatlas-tsv data/peptideatlas/peptideatlas-phospho-202512.tsv

# Build from both sources
hvantk ptm build -o data/ptm/ --output-ht ptm.ht \
    --ptm-tsv data/uniprot/uniprot-ptm.tsv \
    --peptideatlas-tsv data/peptideatlas/peptideatlas-phospho-202512.tsv
```
