# ClinVar Query Examples

Query and filter ClinVar variant annotations using hvantk's `ClinVarVariantTableStreamer`.

## Prerequisites

Build a ClinVar Hail Table first. The ClinVar plugin has a built-in downloader, so this single command downloads the latest VCF into `data/clinvar/` and builds the table in one go:

```bash
hvantk reprocess clinvar:variants \
  --raw-dir data/clinvar/ \
  --output clinvar.ht \
  --plugin-arg reference_genome=GRCh38
```

If you already have `clinvar.vcf.bgz` (or `.vcf.gz`), place it under `data/clinvar/` and add `--skip-download`.

## Basic Usage

The streamer wraps a built ClinVar `AnnotationTable` and exposes typed query
methods. Construct it from a saved table with `from_path`, then drop to the
underlying Hail Table with `to_hail()`:

```python
from hvantk.skills.clinvar.streamers import ClinVarVariantTableStreamer

streamer = ClinVarVariantTableStreamer.from_path("clinvar.ht")

ht = streamer.to_hail()
print(f"Total variants: {ht.count()}")
```

## Filtering by Gene Set

`filter_to_genes` keeps variants whose `info.GENEINFO` gene matches the set and
returns a Hail Table:

```python
brca_ht = streamer.filter_to_genes({"BRCA1", "BRCA2", "TP53"})
print(f"Variants in panel: {brca_ht.count()}")
```

## Filtering by Pathogenicity

`filter_by_pathogenicity` filters on the ClinVar `info.CLNSIG` clinical
significance labels:

```python
pathogenic_ht = streamer.filter_by_pathogenicity(
    ["Pathogenic", "Likely_pathogenic"]
)
print(f"Pathogenic variants: {pathogenic_ht.count()}")
```

## Aggregating Results

The methods return ordinary Hail Tables, so aggregate them natively:

```python
import hail as hl

brca_ht = streamer.filter_to_genes({"BRCA1", "BRCA2", "TP53"})
by_sig = brca_ht.aggregate(hl.agg.counter(brca_ht.info.CLNSIG))
print(by_sig)
```

## Deriving a Training-Label Column

`label_training_set` labels missense SNVs 1 (pathogenic) or 0 (benign) from the primary
ClinVar classification -- the text before the first `|` of `CLNSIG`, so
`Pathogenic|risk_factor` counts as pathogenic -- and adds the review stars derived from
`CLNREVSTAT` (0-4). Variants with any other classification are dropped. The library
function behind it, `label_clinvar_variants`, keeps the unlabelled rows as well and
takes the label sets, the consequence and the SNV filter as arguments.

```python
labeled_ht = streamer.label_training_set()
labeled_ht = labeled_ht.filter(labeled_ht.clinvar_review_stars >= 1)
print(labeled_ht.aggregate(hl.agg.counter(labeled_ht.label)))
```

To build a full feature dataset from these labels, see `hvantk training-set build`.

## Tips

- **Data freshness**: Update ClinVar tables regularly (monthly recommended).
- **Gene names**: ClinVar uses HGNC symbols; verify your gene names match.

## Runnable Scripts

See the [`examples/clinvar/`](https://github.com/bigbio/hvantk/tree/main/examples/clinvar/) directory for a complete runnable example.
