# Data Catalog and Hosting Strategy

Because hosting full HT/MT artifacts is heavy, keep a lightweight registry in-repo and host data externally.

## Catalog goals

- Discoverability: What datasets exist, with which versions, and where
- Integrity: Checksums and sizes to verify downloads
- Reproducibility: Parameters and provenance stored alongside

## Format

JSON (default) or YAML (optional). One file can contain many entries or use one-per-dataset under hvantk/resources/catalog.

Entry fields:
- id: unique dataset id (e.g., clinvar.ht, ucsc.pancreas.mt)
- version: semantic or date version
- type: table|matrixtable
- build: GRCh37|GRCh38
- uri: remote object storage or DOI (s3://, gs://, https://zenodo.org/record/...)
- checksum: sha256 or md5
- size_bytes: integer
- n_partitions: integer (optional)
- schema_hash: stable hash of schema
- params: map of builder params (selectors)
- created: ISO8601 timestamp
- provenance:
  - source: original URLs/DOIs
  - builder: module and function used
  - commit: git commit hash

## Example

- id: clinvar.ht
  version: 2025-05-01
  type: table
  build: GRCh38
  uri: https://example.org/hvantk/artifacts/clinvar/2025-05-01/clinvar.ht
  checksum: 3f785...abc
  size_bytes: 123456789
  n_partitions: 256
  schema_hash: a1b2c3d4
  params: { reference_genome: GRCh38 }
  created: 2025-05-01T12:00:00Z
  provenance:
    source:
      - https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    builder: hvantk.tables.table_builders.create_clinvar_tb
    commit: deadbeef

## Workflow

- Keep manifests small and versioned in Git
- Publish big artifacts to durable storage (S3/GCS/Zenodo/figshare/OSF)
- Tools resolve catalog: prefix path with catalog:<id>@<version> to fetch/locate
- Verify checksum on download; only then import with Hail

## Minimal resolver

- Given catalog id@version, return a local path. If missing, download to ~/.cache/hvantk/<id>/<version>/ and verify checksum
- Keep the resolver as a tiny utility in hvantk.utils (future work)
# Developing hvantk

This guide makes the developer workflow explicit and addresses heterogeneous omics, limited hosting, and pipeline composition.

## Core ideas

- Data products: Immutable Hail Tables/MatrixTables with a clear schema and key. Treat each as a versioned artifact.
- Slice-first builds: Builders accept selectors (genes/regions, tissues/cell types, timepoints, cohorts) to avoid building full tables.
- Caching by parameters: Persist outputs under a path that encodes dataset id, version, and a short hash of the selectors/params.
- Streamers: Stateless transformers that read a Table/MatrixTable and emit another, composed in order.
- Recipes: JSON by default to wire streamers and parameters reproducibly (YAML optional).

## Contracts

### 1) Dataset builder contract
- Inputs: raw files or external URIs; optional selectors; optional reference genome.
- Output: Hail Table/MatrixTable with:
  - stable key: e.g., Table keyed by [locus, alleles] (variants) or gene_id (genes)
  - typed fields: document types and nullability
  - metadata: a dict with {dataset_id, version, source, build, params}
- Behavior: deterministic for the same input+params; idempotent when overwrite=False.

### 2) Streamer contract (Hail)
- setup(): prepare resources; may init Hail
- stream(): yield chunks when first in pipeline
- process_chunk(ht_or_mt): return transformed chunk
- teardown(): free resources, stop Hail if started
Reference: hvantk/data/data_streamer.py and hvantk/annotation/annotation_streamer.py

### 3) Recipe format
- JSON containing inputs, selectors, streamers list with names+params, and outputs (YAML optional). See docs/STREAMERS_AND_RECIPES.md.

## Folder conventions

- hvantk/tables: builders that create HT/MT from raw
- hvantk/commands: CLI entry points
- hvantk/data and hvantk/annotation: streamers and processors
- hvantk/resources: small JSON manifests (YAML optional), not big data
- data/ (user space): actual HT/MT outputs and caches (ignored in VCS)

## Slice-first builder pattern

- Add selectors: genes, regions, tissues, cell_types, timepoints, cohorts, columns
- Implement efficient filtering early (before heavy joins)
- Name outputs like: data/{dataset_id}/{version}/{param_hash}/table.ht
  - param_hash: short stable hash of a canonical JSON of params

## Minimal example: add a new gene constraint dataset

- Builder: hvantk/tables/table_builders.py (follow create_gnomad_constraint_gene_metrics_tb)
- CLI: add a flag in hvantk/commands/make_table_batch_cli.py to call the builder
- Manifest: add an entry to hvantk/resources/catalog.yaml (see DATA_CATALOG.md)
- Tests: add a tiny TSV to hvantk/tests/testdata and test import/select/key

## Testing

- Run: pytest -q
- Add: small synthetic fixtures under hvantk/tests/testdata
- Prefer 3 kinds of tests:
  - schema tests: keys present, field types
  - functional tests: selectors filter correctly
  - pipeline tests: 2-3 streamers compose on tiny data

## CLI patterns

- hvantk ucsc-downloader: fetch raw data
- hvantk ucsc-matrix: build a small MatrixTable from raw
- hvantk mktable-batch: build a bundle of Table artifacts from local raw via JSON recipe (YAML optional)
Extend these commands rather than creating many new ones.

## Reproducibility

- Record parameters: write a sidecar JSON (params + git commit + builder version) next to outputs
- Use immutability: never overwrite by default; produce a new version or param_hash directory
- Keep data sources referenced by DOI/URI + checksum in the catalog

## Performance tips

- Repartition early for large imports
- Select only needed fields before checkpoints
- Use write()/checkpoint() judiciously; avoid too many small partitions
