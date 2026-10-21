# Gene sets

`hvantk genesets` extracts or prepares `GeneSetCollection` JSON files: named sets of
genes with a shared background universe. Three subcommands (`clingen`, `gencc`,
`cosmic`) extract sets from a curated gene-disease table you've already built;
`prepare` converts a plain-text gene panel into the same format.

## When to use it

Use it whenever an analysis tool needs a gene set as input rather than a single
annotation column: [`hvantk psroc`](psroc.md)'s `--gene-sets`,
[`hvantk enrichex`](enrichex.md)'s `burden -s` / `overlap -s`, or a
[cohort manifest](cohort.md)'s `labels.gene_set`.

## Inputs

- `clingen` / `gencc` / `cosmic` each take `--ht`, the provider's built gene-disease
  table (`hvantk reprocess clingen:gene-disease` / `gencc:submissions` /
  `cosmic-cgc:submissions`).
- `prepare` takes `-i`, a headerless two-column TSV: `gene_set_name<TAB>gene_symbol`,
  one membership per line (`#` comments and blank lines ignored). It expects
  **HGNC-style human gene symbols** -- see [Caveats](#caveats).

## Worked example

```bash
# Extract ClinGen gene sets by Gene Curation Expert Panel
hvantk genesets clingen \
  --ht out/clingen_gene_disease.ht \
  --group-by gcep --min-genes 20 \
  -o out/clingen_gcep.json

# Convert a plain panel TSV, validating symbols against HGNC
hvantk genesets prepare \
  -i panels.tsv -o out/panels.json \
  --hgnc out/hgnc_lookup.ht --export-gmt out/panels.gmt
```

## Options reference

All four subcommands share `-o, --output PATH` (`.json` or `.gmt`; `.gmt` loses the
metadata JSON carries) and `--overwrite` (without it, an existing output is a hard
error).

### `genesets clingen` / `gencc` / `cosmic`

| Option | Required | Description |
| --- | --- | --- |
| `--ht PATH` | yes | The provider's built gene-disease table. |
| `--group-by [...]` | no | `clingen`: `gcep` (default), `disease`, `keyword`. `gencc`: `submitter` (default), `disease`, `keyword`. `cosmic`: `tumour-type` (default), `role`, `tissue`, `keyword`. |
| `--min-classification [...]` | no | `clingen`: `Definitive`/`Strong`/`Moderate` (default)/`Limited`/`Disputed`/`Refuted`. `gencc`: `Definitive`/`Strong`/`Moderate` (default)/`Supportive`/`Limited`/`Disputed Evidence`/`Refuted Evidence`. `cosmic`: `Tier 1` (default) / `Tier 2`. |
| `--mutation-context [somatic\|germline\|both]` | `cosmic` only | Default `both`. |
| `--min-genes INTEGER` | no | Drop groups below this size. Default `0`. |
| `--categories-json PATH` | only with `--group-by keyword` | JSON mapping a category name to a list of disease-label keywords (passed to `aggregate_by_disease_category`). |

`--group-by gcep`/`submitter`/`tumour-type`/`role`/`tissue` groups directly on the
table's own field; `disease` groups on each individual disease label; `keyword`
aggregates diseases into the categories you declare.

### `genesets prepare`

| Option | Required | Description |
| --- | --- | --- |
| `-i, --input PATH` | yes | The panel TSV. |
| `--hgnc PATH` | no | A built HGNC Hail Table (`.ht`) for symbol validation and alias resolution. |
| `--min-genes INTEGER` | no | Drop gene sets below this size, after HGNC validation. Default `0`. |
| `--background PATH` | no | Text file, one symbol per line; default is the union of genes across all sets. |
| `--export-gmt PATH` | no | Also write a `.gmt` alongside the JSON. |

## The `GeneSetCollection` JSON format

Verified against `hvantk/core/utils/gene_sets.py` (`GeneSetCollection.to_dict`/`from_dict`):

```json
{
  "gene_sets": {
    "T_cell": {"name": "T_cell", "genes": ["CD3D", "CD3E", "CD3G"], "source": "...", "n_genes": 3, "metadata": {}}
  },
  "background_genes": ["CD19", "CD3D", "CD3E", "CD3G", "CD79A", "CD79B"],
  "n_gene_sets": 1,
  "n_background": 6,
  "source_description": "...",
  "metadata": {}
}
```

`genes` and `background_genes` are always sorted lists, for a reproducible diff. The
`clingen`/`gencc`/`cosmic` extractors and `prepare` all set each `GeneSet.source` and
the collection's `source_description` to the same string (e.g. `clingen_gcep`); the
extractors always default `background_genes` to the union of genes across the sets
*they produced* -- only `prepare --background` lets you supply an external universe.
`prepare` additionally stamps `metadata` with `created_by`, `input_file`,
`hgnc_validated`, and (when `--hgnc` was given) `aliases_resolved` /
`unrecognized_symbols`.

### GMT export

`.gmt` is `<set_name><TAB><description><TAB><gene1><TAB><gene2>...`, one line per set,
sorted by name; `description` is a chosen metadata field or the set's `source`, empty
string if neither is set. GMT has no field for the background universe or collection
metadata -- round-tripping through `.gmt` loses both (`load_gmt`'s background then
defaults to the union of genes in the file, same rule as the JSON path).

## HGNC validation

`prepare --hgnc` classifies every symbol as canonical, a resolvable alias/previous
symbol (rewritten to the canonical form, de-duplicating within its set), or
unrecognized (kept as-is, with a warning -- it will not match any gene downstream).
Without `--hgnc`, `prepare` still hard-rejects a symbol it can tell isn't a human gene
symbol at all before writing anything (see [Caveats](#caveats)); it just skips alias
resolution.

## Consumers

| Consumer | How |
| --- | --- |
| `hvantk psroc` | `--gene-sets <collection.json>` |
| `hvantk enrichex burden` / `overlap` | `-s <collection.json>` |
| `hvantk cohort` | a manifest's `labels.gene_set` (see [External cohorts](cohort.md)) |
| `hvantk expression markers` | *produces* one (see [Expression analysis](expression.md)) |

## Python API

```python
from hvantk.core.utils.geneset_io import parse_geneset_tsv, validate_with_catalog
from hvantk.core.utils.gene_sets import load_gene_sets_from_dict, GeneSetCollection
from hvantk.skills.hgnc.streamers import HGNCGeneCatalogStreamer

result = parse_geneset_tsv("panels.tsv")        # ParseResult(gene_sets, warnings, ...)

catalog = HGNCGeneCatalogStreamer.from_path("out/hgnc_lookup.ht")
vr = validate_with_catalog(result.gene_sets, catalog)  # ValidationResult

collection = load_gene_sets_from_dict(vr.gene_sets, source="prepare-geneset")
collection.save("panels.json")
collection.save_gmt("panels.gmt")

reloaded = GeneSetCollection.load("panels.json")  # or .load_gmt(...)
```

`parse_geneset_tsv` and `validate_with_catalog` live in `hvantk.core.utils.geneset_io`;
everything else -- `GeneSet`, `GeneSetCollection`, `load_gene_sets_from_dict`,
`load_gene_set` (a plain-text single set), `load_gene_sets` (format-sniffing
`.json`/`.gmt` loader), `load_marker_genes`, `extract_marker_gene_sets` -- lives in
`hvantk.core.utils.gene_sets`.

## The ClinGen / GenCC streamer exports

The query/aggregation base class (`hvantk.core.streamers.gene_disease_table.GeneDiseaseTableStreamer`)
backs `ClinGenGeneDiseaseTableStreamer`, `GenCCGeneDiseaseTableStreamer`, and
`CosmicCGCGeneDiseaseTableStreamer` alike, so `export_for_enrichex`, `to_gene_set`, and
ID translation work the same way regardless of provider:

```python
from hvantk.skills.clingen.streamers import ClinGenGeneDiseaseTableStreamer

streamer = ClinGenGeneDiseaseTableStreamer("out/clingen_gene_disease.ht")

definitive = streamer.get_genes_by_classification("Definitive")

cancer_genes = streamer.get_genes_by_disease(
    ["cancer", "carcinoma", "tumor"], match_mode="contains", min_classification="Moderate",
)

stats = streamer.compute_stats()
summary = streamer.classification_summary()

streamer.export_for_enrichex("out/clingen_gene_sets.json", min_classification="Moderate")

# Translate to Ensembl gene IDs using the HGNC catalog
from hvantk.skills.hgnc.streamers import HGNCGeneCatalogStreamer

catalog = HGNCGeneCatalogStreamer.from_path("out/hgnc_lookup.ht")
ensembl_ids = streamer.get_genes_by_classification(
    min_classification="Definitive", gene_catalog=catalog, output_id_type="ensembl_gene_id",
)
hgnc_ids = streamer.to_gene_set(
    min_classification="Moderate", gene_catalog=catalog, output_id_type="hgnc_id",
)
```

GenCC's provider-specific grouping is by submitter instead of GCEP:

```python
from hvantk.skills.gencc.streamers import GenCCGeneDiseaseTableStreamer

streamer = GenCCGeneDiseaseTableStreamer("out/gencc_submissions.ht")
by_submitter = streamer.get_geneset_per_submitter(min_classification="Moderate", min_genes=20)
consensus = streamer.consensus_genes(min_submitters=2)  # genes >=2 submitters agree on
```

`output_id_type` accepts whatever the `gene_catalog` you pass supports; for
`HGNCGeneCatalogStreamer` that is `hgnc_id`, `gene_symbol`, `ensembl_gene_id`,
`entrez_id`, or `uniprot_id`.

## Caveats

- **`hvantk.core.models.GeneSet` (one of the four dataset types) is unrelated to
  `GeneSetCollection`'s `GeneSet`** (`hvantk.core.utils.gene_sets.GeneSet`), despite the
  shared name. Everything on this page is the latter. No in-tree plugin currently
  builds the dataset-type `GeneSet` (`.geneset.json`).
- `prepare` (with or without `--hgnc`) hard-rejects a panel containing Ensembl gene or
  transcript IDs, bare numeric (Entrez-looking) tokens, or whitespace inside a symbol,
  and rejects the whole file if more than half of its unique symbols look like mouse
  gene symbols (`Abc1`-shaped) -- it is a human-HGNC-symbol format, not a general
  identifier container.
- `--group-by keyword` (on any of `clingen`/`gencc`/`cosmic`) requires
  `--categories-json`; without it the command exits before touching the table. On
  `cosmic`, `--mutation-context` applies in every `--group-by` mode, including
  `keyword`.

## See also

- [External cohorts](cohort.md) -- a cohort manifest's `labels.gene_set`
- [Expression analysis](expression.md) -- `markers` produces the same JSON format
- [PS-ROC](psroc.md), [EnrichEx](enrichex.md) -- the analysis tools that consume it
- [Building datasets](../guide/usage.md) -- build the ClinGen/GenCC/COSMIC tables this extracts from
