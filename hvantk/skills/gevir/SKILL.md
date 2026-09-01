# gevir

GeVIR (Gene Variation Intolerance Rank) provides gene-level metrics quantifying
the intolerance of human genes to functional variation. The resource was
published in PMID 31873297 (Abramovs, Brass & Tassabehji, 2020, *Nature
Genetics* 52(1):35-39; DOI 10.1038/s41588-019-0560-2) and ranks **19,361**
protein-coding genes by their intolerance to variation, derived from the density
and spatial distribution of protein-coding variants observed across ~138,632
gnomAD exome and genome sequences. GeVIR is a gene-level metric — it is **not** a
variant-level pathogenicity score.

The metric table is distributed as **Supplementary Table 2** of that paper, served
from Springer's static-content CDN as the article's MOESM3 object (an `.xlsx`
workbook, ~10.3 MB; only that one of the six MOESM slots is public). The authors'
repository at https://github.com/gevirank/gevir ships the **analysis code only** --
its `tables/` directory holds a placeholder file -- so it is not a source for the
table and cannot be used as a drift target.

## Dataset

- `gevir:metrics` — per-gene GeVIR / VIRLoF intolerance rank metrics, keyed by gene_id

The builder is `build_gevir_metrics` in `hvantk/skills/gevir/builder.py`, with
signature `(parsed_input, ctx, **params) -> AnnotationTable`. It imports the
GeVIR TSV inline via `hl.import_table(..., key="gene_id")`, optionally selects a
`fields` subset, and wraps the result with provenance under schema
`gevir-metrics-v1`. The dataset is resolved from `plugin.yaml` by the plugin
loader (`hvantk/core/plugin/loader.py`) and built through `run_builder_for_spec`
(`hvantk/core/plugin/run_builder.py`); there is no separate builder registry.

## Build

```bash
hvantk reprocess gevir:metrics \
  --raw-dir <dir-containing-gevir-tsv> \
  --output <out.ht>
```

Optional plugin args (e.g. field selection) can be passed with
`--plugin-arg fields=...`.

## Notes

The drift probe (`hvantk/skills/gevir/drift_probe.py`, `fetch_fingerprint`) issues
a single HEAD against the article's supplementary object (§ above) and compares
its content-hash ETag plus Content-Length; the workbook body is never
transferred. It replaced a documentation-only stub once that URL was confirmed
addressable (issue #177, which had recorded the source as publication-only and
therefore unprobeable).

**The probe target and the build input are different files.** The probe watches
the upstream `.xlsx` (~10.3 MB); the builder reads a bgzipped TSV
(`gevir_metrics_pmid31873297.tsv.bgz`, ~1-2 MB) derived from sheet `table_2` of
that workbook. The catalog `files` entry describes the derived TSV, which is why
it carries no download URL: fetching the upstream URL does not yield that file
without an extract-and-convert step.

No downloader is implemented yet. A real one would have to extract sheet
`table_2` and BGZF-compress it, not just fetch the URL, so it is more than the
usual thin wrapper — a recommended follow-up. Until then, upstream files are
expected to be materialized externally.

## Schema

The table is keyed by `gene_id` (Ensembl gene ID). The metric is gene-keyed and
build-agnostic; the underlying gnomAD data is v2 (GRCh37).

| field | type | description |
| --- | --- | --- |
| `gnomad_gene_name` | str | HGNC gene symbol as used in gnomAD. |
| `gene_id` | str | Ensembl gene identifier (`ENSG...`); the table key. |
| `canonical_transcript` | str | Ensembl canonical transcript (`ENST...`) used to compute the metric. |
| `gevir_pct` | float64 | GeVIR percentile rank (0-100). Lower percentiles denote genes **more** intolerant to variation (more constrained), derived from the density/spatial clustering of protein-coding variants across gnomAD sequences. |
| `virlof_pct` | float64 | VIRLoF percentile rank (0-100). A combined rank integrating GeVIR with the gnomAD LOEUF loss-of-function constraint metric; lower percentiles denote stronger overall intolerance. |

## Tests

A conformance test in `tests/test_gevir.py` exercises the build via
`run_builder_for_spec` against the bundled fixture at
`hvantk/tests/testdata/raw/gevir/gevir_metrics_pmid31873297.tsv.bgz`.
