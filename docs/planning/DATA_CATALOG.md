# Data Catalog

This document describes the data catalog strategy for hvantk, including versioning, provenance tracking, and dataset hosting.

## Overview

The hvantk data catalog is a lightweight registry that tracks:
- Dataset metadata and provenance
- Version information and checksums
- Remote storage locations (S3, GCS, Zenodo, DOI)
- Processing history and transformations

## Catalog Format

The catalog uses JSON format (YAML optional) to maintain maximum compatibility.

### Catalog Structure

```json
{
  "version": "1.0",
  "datasets": [
    {
      "id": "clinvar_grch38_2024_12",
      "name": "ClinVar Annotations",
      "description": "Clinical significance annotations for genetic variants",
      "version": "2024-12",
      "reference_genome": "GRCh38",
      "format": "hail_table",
      "created": "2024-12-15T10:00:00Z",
      "source": {
        "url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/",
        "provider": "NCBI ClinVar",
        "license": "Public Domain",
        "citation": "doi:10.1093/nar/gkaa1043"
      },
      "files": {
        "hail_table": {
          "uri": "s3://hvantk-data/clinvar/2024-12/clinvar.ht",
          "checksum": "sha256:abc123...",
          "size_bytes": 1234567890
        },
        "raw_vcf": {
          "uri": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20241215.vcf.gz",
          "checksum": "md5:def456...",
          "size_bytes": 987654321
        }
      },
      "schema": {
        "key": ["locus", "alleles"],
        "fields": {
          "clinical_significance": "str",
          "review_status": "str",
          "variation_id": "int32"
        }
      },
      "processing": {
        "builder": "hvantk.builders.clinvar_builder",
        "parameters": {
          "reference_genome": "GRCh38",
          "export_tsv": false
        },
        "timestamp": "2024-12-15T12:00:00Z"
      }
    }
  ]
}
```

## Hosting Strategy

### Limited Hosting Approach

Due to the large size of multiomics datasets, hvantk follows a distributed hosting strategy:

1. **Manifests Only**: Host only catalog manifests and small indices in this repository
2. **Remote Data**: Point to immutable URIs on external platforms
3. **Checksums**: Verify data integrity with checksums
4. **Mirrors**: Support multiple mirrors for availability

### Supported Storage Platforms

- **S3**: Amazon S3 buckets (public or authenticated)
- **GCS**: Google Cloud Storage buckets
- **Zenodo**: Long-term archival with DOIs
- **DOI**: Any dataset with a persistent DOI
- **FTP**: Public FTP servers (for source data)
- **HTTP**: Direct HTTP/HTTPS downloads

## Dataset Registry

### Registry Location

- Main catalog: `docs/registry/catalog.json`
- Dataset-specific metadata: `docs/registry/datasets/`

### Registering a New Dataset

1. **Prepare the dataset**:
   ```bash
   hvantk mktable <type> --raw-input <source> --output-ht <output>
   ```

2. **Compute checksums**:
   ```bash
   # For files
   sha256sum dataset.ht > dataset.ht.sha256
   ```

3. **Upload to remote storage**:
   ```bash
   # Example: Upload to S3
   aws s3 cp dataset.ht s3://bucket/path/ --recursive
   ```

4. **Add to catalog**:
   ```json
   {
     "id": "unique_dataset_id",
     "name": "Dataset Name",
     "version": "1.0",
     "files": {
       "hail_table": {
         "uri": "s3://bucket/path/dataset.ht",
         "checksum": "sha256:..."
       }
     }
   }
   ```

## Versioning

### Version Scheme

Datasets use semantic versioning with date stamps:
- `YYYY-MM` for monthly releases (e.g., ClinVar, gnomAD)
- `v1.0.0` for custom builds
- Include data freeze date in metadata

### Version Management

```json
{
  "id": "dataset_name",
  "versions": [
    {
      "version": "2024-12",
      "release_date": "2024-12-15",
      "uri": "s3://bucket/dataset/2024-12/",
      "deprecated": false
    },
    {
      "version": "2024-11",
      "release_date": "2024-11-15",
      "uri": "s3://bucket/dataset/2024-11/",
      "deprecated": true,
      "deprecation_reason": "Superseded by 2024-12"
    }
  ]
}
```

## Provenance Tracking

Track the full lineage of processed datasets:

```json
{
  "provenance": {
    "source_datasets": [
      {
        "id": "raw_clinvar_2024_12",
        "uri": "https://ftp.ncbi.nlm.nih.gov/...",
        "checksum": "md5:..."
      }
    ],
    "processing_steps": [
      {
        "step": "download",
        "timestamp": "2024-12-15T09:00:00Z",
        "tool": "wget",
        "version": "1.21"
      },
      {
        "step": "build_table",
        "timestamp": "2024-12-15T10:00:00Z",
        "tool": "hvantk",
        "version": "0.1.0",
        "command": "hvantk mktable clinvar ..."
      }
    ],
    "environment": {
      "hail_version": "0.2.132",
      "python_version": "3.10.12",
      "os": "Ubuntu 22.04"
    }
  }
}
```

## Using the Catalog

### Accessing Datasets

```python
from hvantk.catalog import load_dataset, list_datasets

# List available datasets
datasets = list_datasets()

# Load a specific dataset
ht = load_dataset('clinvar_grch38_2024_12')

# Load a specific version
ht = load_dataset('gnomad_constraints', version='v4.1')
```

### Command Line

```bash
# List datasets
hvantk catalog list

# Show dataset info
hvantk catalog info clinvar_grch38_2024_12

# Download dataset
hvantk catalog download clinvar_grch38_2024_12 --output /data/
```

## Caching Strategy

### Parameter-based Caching

Cache processed data by parameter hash to enable reuse:

```python
# Cache key based on parameters
cache_key = hash({
  'dataset': 'expression_atlas',
  'tissue': 'heart',
  'timepoint': 'adult',
  'min_expression': 5.0
})

# Check cache
cached = check_cache(cache_key)
if cached:
    return load_cached(cache_key)
```

### Cache Location

- Local: `~/.hvantk/cache/`
- Shared: `/shared/hvantk/cache/` (for clusters)
- Remote: S3/GCS buckets

## Slice-First Approach

For heterogeneous multiomics data, prefer slicing before building:

```python
# Good: Build only what's needed
build_expression_table(
    genes=['BRCA1', 'TP53'],
    tissues=['heart', 'brain'],
    timepoints=['adult']
)

# Avoid: Building the full dataset
build_expression_table()  # Everything!
```

## Best Practices

1. **Use immutable URIs**: Never modify published datasets
2. **Include checksums**: Always verify data integrity
3. **Version everything**: Track versions for reproducibility
4. **Document sources**: Include citations and licenses
5. **Test before publishing**: Validate datasets before adding to catalog
6. **Provide examples**: Include usage examples for each dataset

## See Also

- [Developer Guide](DEVELOPING.md) - Development workflow
- [Streamers and Recipes](STREAMERS_AND_RECIPES.md) - Data processing pipelines
- [Annotation Sources](../library/annotation-sources.md) - Available data sources
- [Registry](../registry/) - Dataset registry files
