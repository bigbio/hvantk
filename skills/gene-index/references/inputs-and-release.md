# The two inputs, the Ensembl pin, and the receipt

## The two inputs

`annotate spine` takes two already-BUILT Hail tables, never raw files:

- **`ensembl-gene:structure`** (`hvantk/skills/ensembl_gene/structure/`): the canonical
  per-gene table, keyed on `gene_id`, parsed from the pinned-release Ensembl GTF. No
  parse stage is declared -- the builder reads the downloaded `.gtf.gz` directly (a
  pure-pandas parse, then `hl.Table.from_pandas`). It carries `gene_name`, `chromosome`,
  `gene_start`, `gene_end`, `gene_biotype`, `cds_length`, `n_coding_exons`,
  `n_transcripts`, `mane_select` (a plain `str`, the MANE Select transcript id if any,
  empty otherwise) and `cds_transcript`. Acquisition:
  ```bash
  hvantk download ensembl-structure --raw-dir <raw dir> [--overwrite]
  hvantk reprocess ensembl-gene:structure --raw-dir <raw dir> \
    --output <out>.ht --skip-download --no-check-drift
  ```
  The download is a single `.gtf.gz` file, tens of megabytes -- not the much larger
  figure some dataset-size estimates quote for a full Ensembl release.

- **`hgnc:lookup`** (`hvantk/skills/hgnc/`): the HGNC complete-set TSV as a Hail table,
  keyed on `hgnc_id`, carrying (among many columns) `ensembl_gene_id` and `gene_group`
  (`array<str>`, a gene can sit in more than one HGNC group) and its OWN `mane_select`
  (`array<str>` -- a different shape from structure's scalar `mane_select`; see below).
  No parse stage either. Acquisition:
  ```bash
  hvantk download hgnc --output <raw dir>/hgnc_complete_set.txt [--overwrite]
  hvantk reprocess hgnc:lookup --raw-dir <raw dir> \
    --output <out>.ht --skip-download --no-check-drift
  ```
  Unlike the Ensembl downloader, which silently returns the cached path when the target
  already exists, `hvantk download hgnc` raises if its `--output` file exists and
  `--overwrite` is not given -- check for an existing file yourself before calling it
  unconditionally in a script.

Both builds write a sidecar `<output>.provenance.json` next to the table
(`hvantk/core/io/_manifest.py`), with fields `plugin`, `dataset`, `plugin_version`,
`source_fingerprint`, `schema_id`, `build_timestamp`, `builder_commit`, `parents` and
`build_parameters`. `builder_commit` is how you confirm two inputs (or an input and the
checkout building the gene index) come from the same checkout.

## The Ensembl release pin

One constant, `hvantk/resources/ensembl_release.py::ENSEMBL_RELEASE` (currently `"113"`),
fixes the GTF filename (`Homo_sapiens.GRCh38.{RELEASE}.gtf.gz`) and download URL for the
whole toolkit -- not just `ensembl-gene:structure`, but also the PTM coordinate mapper.
It lives in `resources/` rather than in the `ensembl_gene` skill because both a skill and
an algorithm depend on it, and the one-way dependency rule (`skills`/`algorithms` ->
`core`/`resources`) forbids either from depending on the other.

Bumping the pin means: edit this one file, re-download and rebuild
`ensembl-gene:structure` from the new release, rebuild the gene index from that new
build, and then rebuild everything joined onto the gene index (every `annotate prepare`
output, the composed Layer-1 matrix, anything built from it) -- a gene's CDS length,
MANE Select assignment and coordinates can all change release to release, and two
definitions that drift apart is a silent correctness bug, not a loud one. There is no
GRCh37 gene index and no genome argument on `annotate spine`; a GRCh37 cohort must be
lifted over before it can join.

## HGNC has no version pin

HGNC publishes a rolling "complete set" dump with no versioned URL: `hvantk download
hgnc` always fetches whatever `genenames.org` serves today. Re-running the download on a
different day can therefore change the row count and the gene index's HGNC mapping rate
even with no code change on your side. If a mapping rate looks different from a previous
run, check whether the HGNC raw file is actually the same vintage before suspecting a
code or biotype problem.

## Why the gene index itself has no sidecar

`annotate spine`'s command (`hvantk/tools/annotation/annotate_cli.py::spine_cmd`) writes
its output with a bare `ht.write(output, overwrite=True)`, not through
`hvantk.core.io.save` -- the one call that stamps a `.provenance.json` manifest
elsewhere in the toolkit. The gene index is the deliberate exception: it has no sidecar
of its own. `templates/describe-gene-index.py`'s receipt is what replaces it -- a content
digest of the gene index itself (`content_sha256`), plus the sha256 of the two input
datasets' OWN `.provenance.json` sidecars (not their full contents, just a digest of
each), so the receipt answers "what exactly is this gene index, and what two inputs made
it" without needing a sidecar on the gene index to do it.

## `mane_select`: structure's field, not HGNC's

The gene index's `mane_select` column comes from `ensembl-gene:structure`'s own
`mane_select` (a plain `str`), never from HGNC's `mane_select` (`array<str>`).
`build_spine` (`hvantk/algorithms/annotation/spine.py`) only pulls `hgnc_id` and
`gene_group` out of the HGNC table on purpose -- pulling HGNC's own `mane_select` too
would collide with structure's field of the same name but a different, incompatible
type.
