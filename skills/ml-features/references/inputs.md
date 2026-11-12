# Inputs of a training set: staging and building them

`hvantk training-set build` reads only built datasets. Every `--input NAME=PATH` it takes
is a Hail table built by `hvantk reprocess`, the gene matrix is an `annotate compose`
output, and a `layer1` axis that carries a `matrix:` block is prepared from an expression
summary (`.h5ad`). This page lists, for the shipped example spec
(`examples/training_sets/clinvar_heart.yaml`), where each input comes from, how its raw
file is staged and how it is built. A spec of your own needs the subset it names: list
them with `grep -nE '^\s+(- )?(input|hgnc_input|source):' <spec.yaml>`.

Every command below runs as the skill's Procedure says: downloads on the login node
through `skills/hpc-hail-job/templates/run-in-image.sh`, builds as jobs through
`skills/hpc-hail-job/templates/submit-hail-job.sh`, each starting with
`source "$HVANTK_SITE_PROFILE" &&`. `RAW` and `OUT` stand for absolute directories under
`$HVANTK_WORK`.

## Variant-level inputs (the spec's `label`, `variant` and `gene` sections)

| `--input` name (example) | Dataset | Stage the raw file | Build |
| --- | --- | --- | --- |
| `clinvar` | `clinvar:variants` | `hvantk download clinvar --version YYYYMMDD --output-dir RAW/clinvar` | `hvantk reprocess clinvar:variants --raw-dir RAW/clinvar --output OUT/clinvar.ht --skip-download --no-check-drift` |
| `dbnsfp` | `dbnsfp:variants` | bring your own (below) | `hvantk reprocess dbnsfp:variants --raw-dir RAW/dbnsfp --output OUT/dbnsfp.ht --skip-download --no-check-drift` |
| `insider_variants` | `insider:variants` | manual: the BED's URL is in `hvantk catalog show INSIDER_v1.0` | `hvantk reprocess insider:variants --raw-dir RAW/insider --output OUT/insider_variants.ht --skip-download --no-check-drift` |
| `ptm` | `uniprot-ptm:sites` | `hvantk download uniprot-ptm --output-dir RAW/uniprot_ptm`, plus the pinned Ensembl GTF (below) | `hvantk reprocess uniprot-ptm:sites --raw-dir RAW/uniprot_ptm --output OUT/ptm_sites.ht --intermediate OUT/ptm_sites.intermediate.tsv.bgz --skip-download --no-check-drift` |
| `hgnc` | `hgnc:lookup` | built by the `gene-index` skill | the same `hgnc-lookup.ht` the gene index was built with |

**ClinVar: always pin `--version`.** `latest` changes every week, and a rebuild from a
different release is a different dataset (other rows, other labels). Use the dated
release (`YYYYMMDD`, a weekly archive name) and keep it in your run notes. The sidecar's
`source_fingerprint` and the manifest's `input_provenance` record which file was read.

**dbNSFP is bring-your-own.** It has no downloader: the academic release is obtained by
registration with the dbNSFP maintainers (the plugin's resource spec, `hvantk/skills/dbnsfp/SKILL.md`,
section 2, has the route). The raw directory must hold exactly one `.bgz` file (the
builder picks the one `*.bgz` it finds). The built table of the full release is about
70 GB; build it once and reuse it.

The release ships one plain-gzip file per chromosome. `bgzip` may be on neither the
login node nor the image, so make the one `.bgz` without it: concatenate the files on the
login node (`cat dbNSFP4.9a_variant.chr*.gz > dbNSFP4.9a_variant.gz`; gzip files
concatenated are one gzip stream), then convert it as a job, `hvantk utils convert-bgz
dbNSFP4.9a_variant.gz -o RAW/dbnsfp/dbNSFP4.9a_variant.bgz` (single-threaded without
`bgzip`; allow several hours). The header line then repeats once per chromosome, which
the builder skips.

**dbNSFP rankscores must be numbers.** The current builder parses every `*_rankscore`
field to `float64`. A table built by an older checkout keeps them as strings, and
`training-set build` passes a string column through without any error: the dataset then
carries text where a model expects numbers. Before using a dbNSFP table built
elsewhere, check one field's type and rebuild if it is `str`:

```bash
source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
  --name <run>-dbnsfp-types --cpus 2 --mem-gb 8 --time 0-00:10:00 --wait -- \
  python -c "import hail as hl; from hvantk.core.utils.hail_context import init_hail; init_hail(quiet=True); t = hl.read_table('/abs/path/dbnsfp.ht'); print({f: str(t[f].dtype) for f in t.row if f.endswith('_rankscore')})"
```

**INSIDER variants** come from a BED of about 1.2 GB with no downloader (manual
acquisition; `hvantk/skills/insider/SKILL.md`). Keep the file as downloaded: the builder
reads its `track name=` lines to name each interface.

**UniProt PTM** sites are mapped onto the pinned Ensembl GTF
(`hvantk/resources/ensembl_release.py`), which the build downloads into the raw directory
when it is not there. On a site whose compute nodes have no internet
(`HVANTK_EGRESS=blocked`), copy that file into `RAW/uniprot_ptm` before the build job;
the `gene-index` skill's `ensembl-gene` raw directory has it. UniProt keeps no archived
releases, so the download is today's snapshot (`--version` only labels the file): note
its date in the run notes.

## Gene-level inputs (the spec's `layer1` axes)

Each axis is prepared from one input with `annotate prepare --axis <axis>`, and all of them
are composed into the gene matrix (the `gene-index` skill's
`skills/gene-index/references/annotation-axes.md` has both commands).

| Axis (example) | Prepared from | Stage and build |
| --- | --- | --- |
| `gnomad_constraint` | `gnomad-metrics:metrics` (`.ht`) | `hvantk download gnomad-metrics --output RAW/gnomad_metrics/<file>` (default release `v2.1.1`), then `hvantk reprocess gnomad-metrics:metrics --raw-dir RAW/gnomad_metrics ...` |
| `gevir` | `gevir:metrics` (`.ht`) | `hvantk download gevir --output RAW/gevir/<file>`, then `hvantk reprocess gevir:metrics --raw-dir RAW/gevir ...` |
| `insider_interfaces` | `insider:interfaces` (`.ht`) | `hvantk download insider-interfaces --output RAW/insider_interfaces/H_sapiens_interfacesALL.txt`, then `hvantk reprocess insider:interfaces --raw-dir RAW/insider_interfaces ...`; prepare it with `--hgnc` (it is keyed by `uniprot_id`) |
| `heart_bulk` | an expression summary (`.h5ad`) | E-MTAB-6814: `hvantk download expression-atlas --accession E-MTAB-6814 --download-path RAW/emtab6814`, then `hvantk reprocess expression-atlas:dataset --raw-dir RAW/emtab6814 --output OUT/e-mtab-6814.h5ad --skip-download --no-check-drift --plugin-arg level=gene`, then `hvantk expression summarize -m OUT/e-mtab-6814.h5ad --group-by developmental_stage --filter-by organism_part=heart --min-cells 1 -o OUT/heart_bulk_stage_summary.h5ad` |
| `farah` | an expression summary (`.h5ad`) | the Farah 2024 fetal-heart atlas from the UCSC Cell Browser, collection `hoc`, whole heart: `hvantk download ucsc --dataset hoc/all-heart --output-dir RAW/farah` (into `RAW/farah/hoc/all-heart/`; `--list-datasets --search hoc` shows the collection), then `hvantk expression summarize-ucsc -e RAW/farah/hoc/all-heart/exprMatrix.tsv.gz -m RAW/farah/hoc/all-heart/meta.tsv --group-by Region --group-by Age --group-by major_cell_class --min-cells 10 -o OUT/farah_region_age_class_summary.h5ad` |
| `hca` | an expression summary (`.h5ad`) | the Heart Cell Atlas global object (`.h5ad`, downloaded by hand), then `hvantk expression summarize -m RAW/hca/<global>.h5ad --group-by cell_type -o OUT/hca_cell_type_summary.h5ad`; check the group field first with `hvantk expression describe` |

A `matrix:` axis needs the summary's group labels to produce the spec's column names
(`<atlas>_<group>_<stat>`, e.g. `heart_bulk_neonate_mean`); a summary grouped on another
field gives other names, and `compose` or `training-set build` stops with
`gene matrix lacks column(s): ...`.

## Was each input built by this checkout?

Every built input carries a `.provenance.json` sidecar. Its `builder_commit` names the
commit of the checkout a build ran from (when it ran through the `hpc-hail-job`
templates). Compare it with `git -C "$HVANTK_REPO" rev-parse HEAD` for each input:

```bash
for f in OUT/*.ht.provenance.json; do python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print(sys.argv[1], d.get('dataset'), d.get('builder_commit'))" "$f"; done
```

A different commit is not wrong by itself, but a builder change between the two commits
can change the table (the dbNSFP rankscore types are the known case). After the build, the
manifest's `input_provenance` keeps the same fields for every input, so the dataset
records what it was built from.
