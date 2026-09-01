[![Python Package using Conda](https://github.com/bigbio/hvantk/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/bigbio/hvantk/actions/workflows/python-package-conda.yml)
[![Python application](https://github.com/bigbio/hvantk/actions/workflows/python-app.yml/badge.svg)](https://github.com/bigbio/hvantk/actions/workflows/python-app.yml)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.10-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Docs](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://bigbio.github.io/hvantk)

# hvantk

**Hail-based toolkit for multiomics variant annotation and analysis.**

`hvantk` is a modular toolkit that uses [Hail](https://hail.is/) to annotate and analyze variants, genes, proteins, and expression data from heterogeneous omics sources. The library enables multiomics integration to improve the interpretation of genetic variants.

## Installation

```bash
git clone https://github.com/bigbio/hvantk
cd hvantk
poetry install
eval "$(poetry env activate)"
```

**Prerequisites**: Python >=3.10, Hail

### Optional extras

The base install is intentionally lean — plotting, machine-learning, and a few
provider-specific dependencies are opt-in Poetry extras:

```bash
poetry install --extras "ancestry psroc"   # one or more
poetry install --all-extras                # everything
```

For the full table — what each extra pulls in and which commands need it — see
[Installation → Optional features](docs_site/getting-started/installation.md#optional-features-extras).

Verify it works:

```bash
hvantk utils check-install
hvantk --help
```

## Toolkit

| Tool | Description | Command | Docs |
|------|-------------|---------|------|
| **Downloads** | Acquire external datasets (ClinVar, ClinGen, HGNC, etc.) | `hvantk download <source>` | [Data Sources](docs_site/guide/data-sources.md) |
| **Dataset builds** | Build any plugin dataset (download → parse → build → drift check) | `hvantk reprocess <plugin>:<dataset>` | [Usage Guide](docs_site/guide/usage.md) |
| **HGC** | Joint genotyping pipeline (GVCF combining, QC, format conversion) | `hvantk hgc` | [HGC](docs_site/tools/hgc.md) |
| **Ancestry** | Ancestry inference (PCA + Random Forest classification) | `hvantk ancestry-inference` | [Ancestry](docs_site/tools/ancestry.md) |
| **QTL Cascade** | Molecular QTL integration (eQTL + pQTL cascade, colocalization ABF) | `hvantk qtlcascade` | [QTL Cascade](docs_site/tools/qtlcascade.md) |
| **EnrichEx** | Gene set enrichment (overlap testing + rare variant burden) | `hvantk enrichex` | [EnrichEx](docs_site/tools/enrichex.md) |
| **PS-ROC** | Pathogenicity score ROC evaluation against ClinVar labels | `hvantk psroc` | [PS-ROC](docs_site/tools/psroc.md) |
| **PTM** | Post-translational modification variant classification | `hvantk ptm` | [PTM](docs_site/tools/ptm.md) |
| **Expression** | Expression analysis (summarize, marker extraction) | `hvantk expression` | [Usage Guide](docs_site/guide/usage.md) |
| **Annotate** | Build the gene spine, map sources onto it, compose a gene × feature matrix | `hvantk annotate` | [Usage Guide](docs_site/guide/usage.md#gene-level-annotation-matrices) |
| **Cohort** | Validate and attach external cohorts | `hvantk cohort` | [Usage Guide](docs_site/guide/usage.md#external-cohorts) |
| **Re-rank** | Re-rank genes by multi-omic credibility from a YAML config | `hvantk rerank` | [Example](examples/rerank/README.md) |
| **Gene sets** | Extract or prepare gene set collections (ClinGen, GenCC, COSMIC) | `hvantk genesets` | [Usage Guide](docs_site/guide/usage.md) |

Plus the registry and operational commands: `hvantk plugins` and `hvantk tools` (inspect the
plugin and tool registries), `hvantk catalog` (search the aggregated dataset catalog),
`hvantk drift` (compare a plugin's live source fingerprint against the committed one),
and `hvantk utils` (format conversion, BGZF validation, install diagnostics).

## Architecture

`hvantk` is organized in four code layers (`core/`, `algorithms/`,
`skills/`, `tools/`) plus a substrate-level data registry (`resources/`).
A strict one-way dependency rule is enforced by tests in
[`hvantk/tests/test_dependency_directions.py`](hvantk/tests/test_dependency_directions.py):

<p align="center">
  <img src="docs_site/images/hvantk-platform-architecture.svg" alt="hvantk platform architecture: tools/ depends on skills/, algorithms/, and core/; skills/ depends on algorithms/ and core/; algorithms/ depends on core/. resources/ sits at the substrate level alongside core/ and is consumed by skills/ and tools/. Arrows flow downward only." width="860">
</p>

**Why the directions matter:** `skills/` adapters can rot when upstream APIs
change without algorithms breaking; `algorithms/` evolve without churning
the source-adapter layer. `core/` and `resources/` are the stable substrate
everyone depends on — neither imports upward.

### Data model

Four semantic artifact types live in [`hvantk/core/models/`](hvantk/core/models/),
each backed by one of several native engines:

| Artifact | Backends | On-disk format | Used for |
|---|---|---|---|
| [`AnnotationTable`](hvantk/core/models/annotation_table.py) | `hail` / `pandas` | `.ht/` or `.parquet` | variants, gene-disease pairs, eQTLs, PTM sites |
| [`ExpressionMatrix`](hvantk/core/models/expression_matrix.py) | `anndata` | `.h5ad` | bulk + single-cell expression, proteomics matrices |
| [`VariantMatrix`](hvantk/core/models/variant_matrix.py) | `hail-mt` | `.mt/` | multi-sample variant cohorts (genotypes × samples × multi-field entries) |
| [`GeneSet`](hvantk/core/models/gene_set.py) | (in-memory `frozenset`) | `.geneset.json` | curated gene collections (CHD, MSigDB, …) |

Every artifact carries a [`Provenance`](hvantk/core/models/provenance.py)
record — plugin name, version, source fingerprint, schema id, build
timestamp, and a `parents: tuple[Provenance, ...]` chain for algorithm
derivations. The `@algorithm` decorator stamps input provenances onto
output artifacts automatically, so the build graph is preserved end-to-end.

Artifacts expose a **portable query API** (`filter`, `select`, `join`,
`with_columns`, `group_by().agg()`) via the [`col(...)`](hvantk/core/models/_expr.py)
expression DSL, compiled to either Hail or pandas at execution time —
algorithms can be written backend-agnostically. When the algorithm legitimately
needs the raw native object, [`core_io.load_native(path)`](hvantk/core/io/__init__.py)
returns `(native_obj, Provenance)` zero-cost.

### Plugin contract — adding a new data source

Each data source ships as a self-contained plugin under `hvantk/skills/<plugin>/`,
declared by a [`plugin.yaml`](hvantk/skills/clinvar/plugin.yaml) manifest naming its
builder, drift probe, and downloader. The platform orchestrator
[`run_builder_for_spec`](hvantk/core/plugin/run_builder.py) resolves the manifest,
computes the source fingerprint, calls the builder, validates the returned artifact
against the manifest's `artifact_type` and `schema_id`, and saves it alongside a
sidecar `.provenance.json`. The loader discovers manifests on its own — there is no
registry to edit.

The full contract and the annotated directory tree live in the architecture guide:

- [Plugin contract](docs_site/architecture.md#3-plugin-contract--adding-a-data-source)
  — build sequence diagram, annotated `plugin.yaml`, two-pass loader, streamer
  placement rule
- [Project structure](docs_site/architecture.md#project-structure) — what lives in
  each package, layer by layer

### How to extend

| Add | Where | Pattern |
|---|---|---|
| A new data source | `hvantk/skills/<plugin>/` | Write `plugin.yaml` + `builder.py` (returns one of the four artifact types) + `drift_probe.py`. Loader auto-discovers. |
| A new algorithm | `hvantk/algorithms/<domain>/` | Decorate with `@algorithm(name=…, backends=[…], inputs={…}, outputs={…})`. Operate on Artifact inputs (or use `load_native` for Hail-heavy work). |
| A new CLI command | `hvantk/tools/<domain>/` | Add the click command + a `.tool.yaml` manifest, then an entry in `_LAZY_COMMANDS` in `hvantk/hvantk.py`. Do **not** add a module-level import + `cli.add_command` — that works, but costs every invocation your command's imports (see [Architecture](docs_site/architecture.md#adding-a-new-cli-command)). |
| A new artifact format | `hvantk/core/io/_formats.py` + dispatch in `__init__.py` | Add `save_<artifact>_<ext>` / `load_<artifact>_<ext>`. |

## Documentation

**Full docs site:** [https://bigbio.github.io/hvantk](https://bigbio.github.io/hvantk)

- [Data Sources](docs_site/guide/data-sources.md) -- Available annotations and how to acquire them
- [Examples](docs_site/examples/index.md) -- Tutorials and walkthroughs for each tool
- [Architecture](docs_site/architecture.md) -- Design patterns and extension points

## Citation

If you use hvantk in your research, please cite:

```bibtex
@software{hvantk2024,
  title = {hvantk: Hail-based toolkit for multi-omics variant annotation and analysis},
  author = {Perez-Riverol, Yasset and Audain, Enrique},
  year = {2024},
  url = {https://github.com/bigbio/hvantk}
}
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development workflow, code style, and testing requirements.

```bash
poetry install
pytest -q
hvantk --help
```

## License

MIT License - see [LICENSE](LICENSE).

## Support

- [GitHub Issues](https://github.com/bigbio/hvantk/issues)
- [Documentation](https://bigbio.github.io/hvantk)
