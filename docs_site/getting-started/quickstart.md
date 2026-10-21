# Quick start

This walk-through takes about ten minutes on a laptop. It builds two small public
datasets, turns one of them into gene sets, tests a gene list for enrichment, and
checks a source for upstream changes. You need a working
[installation](installation.md) and internet access. Hail starts a local Spark session
the first time a command needs it, so expect Spark log lines on stderr.

## 1. Check the installation

```bash
hvantk utils check-install
```

It starts Hail, builds a small synthetic MatrixTable and prints its schema. If it fails,
see [Installation → Troubleshooting](installation.md#troubleshooting).

## 2. See what you can build

```bash
hvantk plugins list
```

Every data source is a plugin, and every dataset has an id of the form
`<provider>:<dataset>`, such as `hgnc:lookup` or `clinvar:variants`.
[Data sources](../guide/data-sources.md) says how each one is acquired.

## 3. Build two datasets

```bash
hvantk reprocess hgnc:lookup --raw-dir data/hgnc --output data/hgnc.ht
hvantk reprocess clingen:gene-disease --raw-dir data/clingen --output data/clingen.ht
```

A second run into the same `--raw-dir` stops with `FileExistsError` when the raw file is
already there. Add `--skip-download` to rebuild from it, or use a fresh `--raw-dir`.

Each command downloads the raw file into `--raw-dir`, builds a Hail table at `--output`,
and writes a provenance record next to it (`data/hgnc.ht.provenance.json`):

```json
{
  "plugin": "hgnc",
  "dataset": "hgnc:lookup",
  "plugin_version": "0.2.0",
  "source_fingerprint": "sha256:6f8c05…",
  "schema_id": "hgnc-lookup-v1",
  "build_timestamp": "2026-10-19T22:30:58+00:00",
  "builder_commit": "3af91c0…",
  "parents": [],
  "build_parameters": {}
}
```

Each build ends by probing the upstream source. `drift: clean` means the source still
matches the fingerprint shipped with your hvantk release. `drift: drifted` means it has
changed since (a new release or a larger file); the build still used the current file
and succeeded. See [Check a source for drift](../guide/usage.md#check-a-source-for-drift).

## 4. Turn ClinGen into gene sets

```bash
hvantk genesets clingen --ht data/clingen.ht --min-classification Definitive \
  -o results/clingen_sets.json
```

This writes one gene set per ClinGen expert panel, keeping Definitive gene-disease
pairs only (56 sets in October 2026). The JSON also stores a background universe, by
default the union of all sets in the file (2,101 genes here).
[Gene sets](../tools/genesets.md) covers the other sources and the file format.

## 5. Test a gene list for enrichment

```bash
printf "MYH7\nMYBPC3\nTNNT2\nTNNI3\nTPM1\nACTC1\nMYL2\nMYL3\nLMNA\nSCN5A\n" > my_genes.txt
hvantk enrichex overlap -g my_genes.txt -s results/clingen_sets.json -o results/overlap.tsv
```

`overlap` runs a two-sided Fisher's exact test per gene set against the file's
background and corrects for multiple testing with Benjamini-Hochberg. The ten sarcomere and cardiac genes
come out where you would expect:

```text
gene_set_name                      n_overlap  p_value   p_adjusted  overlap_genes
Hereditary Cardiovascular Disease  5          1.9e-08   1.1e-06     MYBPC3,MYH7,MYL2,SCN5A,TPM1
Dilated Cardiomyopathy             4          1.3e-07   3.5e-06     LMNA,MYH7,SCN5A,TNNT2
Hypertrophic Cardiomyopathy        3          4.3e-06   8.0e-05     ACTC1,MYL3,TNNI3
```

The TSV also carries the odds ratio and its confidence interval.
[EnrichEx](../tools/enrichex.md) explains the test and the case-control burden mode.

## 6. Check a source for drift later

```bash
hvantk drift hgnc:lookup
```

It prints `clean`, or `drifted` with a diff of the fingerprint, and exits 1 when the
source has drifted, so a scheduled job can act on it.

## Where to go next

| You have | Read |
|---|---|
| Per-sample gVCFs to joint-call | [HGC](../tools/hgc.md) |
| A cohort MatrixTable and a reference panel | [Ancestry inference](../tools/ancestry.md) |
| A cohort MatrixTable with cases and controls | [EnrichEx burden](../tools/enrichex.md) |
| A list of genes or variants to annotate | [Building datasets](../guide/usage.md), [Gene annotation matrix](../tools/annotate.md) |
| eQTL, pQTL or GWAS summary statistics | [QTL cascade](../tools/qtlcascade.md) |
| Missense variants near PTM sites | [PTM](../tools/ptm.md) |
| Pathogenicity scores to benchmark | [PS-ROC](../tools/psroc.md) |
| Gene-level cohort results to re-rank | [External cohorts](../tools/cohort.md), [Re-rank](../tools/rerank.md) |
| Expression matrices (`.h5ad`) | [Expression analysis](../tools/expression.md) |
| A Slurm cluster | [Running on HPC](../guide/hpc-migration.md) |
