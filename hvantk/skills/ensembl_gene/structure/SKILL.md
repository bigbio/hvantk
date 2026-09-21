---
name: hvantk:resource-ensembl-gene-structure
description: Per-gene Ensembl structural table (location, biotype, CDS length, MANE Select, transcript count) parsed from the single pinned-release GTF -- Hail Table keyed by gene_id.
status: provisional
backend: hail
domain: mapping
---

# ensembl-gene:structure

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

Provisional. Covers the `ensembl-gene:structure` dataset: the canonical per-gene Ensembl
table -- gene ID, gene name, biotype, coordinates, CDS length, coding-exon count,
transcript count, and MANE Select -- parsed from the single pinned-release Ensembl GTF
(`hvantk/resources/ensembl_release.py`, currently release `113`). Per `builder.py`'s
module docstring, this is the gene spine the annotation pipeline builds on: CDS length
normalises PTM site density, and gene length is a mechanical confounder of rare-variant
burden counts, so both must be available as covariates rather than omitted.

Out of scope for this skill:
- The PTM coordinate mapper's separate GTF parser,
  `hvantk.algorithms.ptm.mapper.parse_ensembl_gtf` -- it keys MANE by gene *symbol* and
  builds per-transcript exon lists for residue-to-genomic coordinate mapping, which this
  dataset's parser (`structure/parse.py`) deliberately does not do; it keys by `gene_id`
  and emits one summary row per gene (see `parse.py`'s module docstring).
- Bumping the Ensembl release pin itself -- an upgrade decision, not something this skill
  automates (see § 8).
- The `.gff3.gz` file listed alongside the GTF in the catalog entry (§ 2) -- this dataset
  reads only the `.gtf.gz` file.

## 2. Source identity

- **Provider:** Ensembl (<https://www.ensembl.org/info/data/ftp/index.html>), genome build
  GRCh38.
- **Catalog entry:** `hvantk/skills/ensembl_gene/catalog/datasets.json`, accession
  `Ensembl_v113` (`data_source: "Ensembl"`, title "Ensembl Gene Annotations", license
  Apache 2.0, `update_frequency: quarterly`). The catalog lists two files
  (`Homo_sapiens.GRCh38.113.gtf.gz` and `Homo_sapiens.GRCh38.113.gff3.gz`); only the
  `.gtf.gz` entry is consumed here.
- **Release pin:** `hvantk/resources/ensembl_release.py` -- `ENSEMBL_RELEASE = "113"`,
  with `ENSEMBL_GTF_FILENAME` and `ENSEMBL_GTF_URL` derived from it as f-strings. Every
  build reads exactly this one release, never "latest". The pin lives in `resources/`
  (substrate, alongside `core`) rather than in this plugin because both this skill and
  `hvantk.algorithms.ptm` consume it, and the dependency rule (`tools -> skills ->
  algorithms -> core`) forbids `algorithms` depending on `skills` (see
  `ensembl_release.py`'s module docstring, and `builder.py`/`parse.py`).
- **Pin consistency is machine-checked in two places**, not just declared:
  - `hvantk/skills/ensembl_gene/tests/test_release_pin.py` (plugin-local) --
    `test_gtf_url_and_filename_carry_the_pinned_release` asserts the release string is
    embedded in the URL/filename constants; `test_catalog_declares_the_same_release`
    asserts the catalog entry's `accession` is `Ensembl_v113` and its `.gtf.gz` file path
    matches `ENSEMBL_GTF_FILENAME`.
  - `hvantk/tests/test_ensembl_release_pin.py` (repo-level) --
    `test_ptm_constants_reuse_the_same_pin` asserts `hvantk.algorithms.ptm.constants`
    re-exports the identical `ENSEMBL_RELEASE` / `ENSEMBL_GTF_URL` /
    `ENSEMBL_GTF_FILENAME`. This lives at the top level (not under
    `hvantk/skills/ensembl_gene/tests/`) specifically because it imports
    `hvantk.algorithms.ptm`, and a file under `skills/` may not import from `algorithms/`
    (`test_dependency_directions`).

## 3. Backend choice + reasoning

`backend: hail`, `domain: mapping` (`plugin.yaml`). The parser itself
(`structure/parse.py`) is deliberately pure Python + pandas -- "no Hail -- so it stays in
the fast test suite", per its module docstring -- but the builder converts the resulting
`pandas.DataFrame` to a Hail Table via `hl.Table.from_pandas(df, key=["gene_id"])`
(`structure/builder.py`). Hail is chosen at that boundary, not for the parse itself,
because downstream consumers (the PTM pipeline, the gene spine) join this table against
other Hail Tables -- the same reasoning `hgnc/SKILL.md` § 3 gives for a small lookup table.

## 4. Raw format & gotchas

Input is a standard Ensembl GTF: tab-separated, 9 columns
(seqname/source/feature/start/end/score/strand/frame/attributes), gzipped as downloaded
(`.gtf.gz`) or plain in test fixtures. The fixture is
`structure/tests/testdata/raw/ensembl-structure/mini.gtf` (14 lines, 3 genes: 2
`protein_coding`, 1 `lncRNA`).

- **Opener:** `gzip.open` if `gtf_path.endswith(".gz")`, else plain `open()`
  (`parse.py`).
- **Comment/header lines:** any empty line or line starting with `#` is skipped -- the
  fixture's first line is `#!genome-build GRCh38` (`parse.py`:
  `if not line or line.startswith("#"): continue`).
- **Feature filtering:** only `gene`, `transcript`, and `CDS` rows are kept; every other
  feature (`exon`, `five_prime_utr`, `start_codon`, ...) is skipped (`parse.py`).
- **Chromosome naming:** Ensembl's bare convention (`1`, `2`, `X`, ... -- no `chr`
  prefix), as seen in the fixture. The parser does **not** call `contig_recoding()`
  (`hvantk/core/utils/genome.py`) or otherwise recode contigs -- `fields[0]` (the seqname
  column) is passed straight through to `chromosome`.
- **Attribute-column parsing is regex-based, not a general key=value parser:**
  `gene_id`, `transcript_id`, `gene_biotype`, and `gene_name` are each pulled from the raw
  attribute string by a dedicated compiled regex (`_RE_GENE_ID`, `_RE_TX_ID`,
  `_RE_BIOTYPE`, `_RE_GENE_NAME` in `parse.py`), each matching `<key> "<value>"`. MANE
  Select is detected by a plain substring check on `transcript` rows,
  `'tag "MANE_Select"' in attrs`, not a tag-list parse.
- **Version-suffix stripping:** both `gene_id` and `transcript_id` have any `.N` version
  suffix removed via `.split(".")[0]` (e.g. `ENSG00000000009.7` ->
  `ENSG00000000009`) -- verified in
  `structure/tests/test_parse.py::test_version_suffixes_are_stripped`.
- **Representative coding transcript:** MANE Select if the gene has one *and* that
  transcript is itself coding (has CDS); otherwise the transcript with the longest summed
  CDS. Ties are broken deterministically on `(cds_length, transcript_id)` tuple
  comparison rather than a bare `max()` over a Python `set` -- the latter would depend on
  the per-process string-hash seed (`parse.py` inline comment). This is pinned by
  `structure/tests/test_parse.py::test_tie_break_deterministic_on_equal_cds_length`, which
  asserts an identical winner across 8 subprocess runs with different `PYTHONHASHSEED`
  values.
- **CDS length:** summed over the representative transcript's `CDS` rows using 1-based
  inclusive GTF coordinates, `end - start + 1` per row.
- **`reprocess` hands the builder a directory, not a file:** this dataset declares
  `lifecycle.download` but no `lifecycle.parse`, so `hvantk reprocess` passes the raw
  directory straight to the builder. `_resolve_gtf_path` (`structure/builder.py`)
  resolves a directory input to `<dir>/<ENSEMBL_GTF_FILENAME>` (the name the downloader
  writes) and returns a file-path input unchanged. This matters because
  `parse_gtf_structure` opens its input with plain `open()`/`gzip.open()`, which raises
  `IsADirectoryError` on a directory -- unlike `hl.import_table`, which tolerates one
  (`builder.py`'s `_resolve_gtf_path` docstring).

## 5. Output contract

Builder returns an `AnnotationTable` (`hvantk.core.models.AnnotationTable`) wrapping a
Hail Table keyed by `gene_id` (string, version suffix stripped), stamped with
`schema_id="ensembl-gene-structure-v1"` via `ctx.provenance(...)`.

**Key:** `gene_id` (Ensembl, version suffix stripped)

| Column | Type | Meaning |
|---|---|---|
| `gene_id` | str | Ensembl gene ID, unversioned |
| `gene_name` | str | Gene symbol read from the `gene` row's `gene_name` attribute |
| `chromosome` | str | Sequence name, Ensembl bare form (e.g. `1`, `X`) |
| `gene_start` | int32 | Gene start (1-based, from the `gene` row) |
| `gene_end` | int32 | Gene end (1-based, inclusive) |
| `gene_biotype` | str | e.g. `protein_coding`, `lncRNA` |
| `mane_select` | str | MANE Select transcript ID, `""` if the gene has none |
| `cds_transcript` | str | Representative coding transcript: MANE if present and coding, else longest-CDS |
| `cds_length` | int32 | Summed CDS bp of `cds_transcript`; `0` for non-coding genes |
| `n_coding_exons` | int32 | CDS feature-row count of `cds_transcript` |
| `n_transcripts` | int32 | Distinct transcript IDs annotated for the gene |

Types and field set are exactly as committed in
`structure/tests/snapshots/schema.json`. Sample rows for the 3-gene fixture (one
MANE-driven pick, one longest-CDS pick, one non-coding gene) are committed in
`structure/tests/snapshots/sample_rows.json`.

**Params:** `protein_coding_only` (bool, default `False`) -- filters the parsed
DataFrame to `gene_biotype == "protein_coding"` rows before the Hail Table is built
(`builder.py`).

## 6. hvantk integration points

- **Plugin manifest:** `hvantk/skills/ensembl_gene/plugin.yaml`, dataset key
  `ensembl-gene:structure`; the loader auto-resolves it via
  `get_registry().get_dataset("ensembl-gene:structure")`.
- **Builder:** `build_ensembl_gene_structure(parsed_input, ctx, **params)` in
  `hvantk/skills/ensembl_gene/structure/builder.py`.
- **Parser:** `parse_gtf_structure(gtf_path) -> pandas.DataFrame` in
  `hvantk/skills/ensembl_gene/structure/parse.py`.
- **Downloader:** `download_dataset(raw_dir, **params)` / Click command `download_cmd`
  in `hvantk/skills/ensembl_gene/structure/cli.py`, wired via `plugin.yaml`'s
  `lifecycle.download` and the `cli:` block (command `ensembl-structure-download`).
  Fetches `ENSEMBL_GTF_URL` into `raw_dir/ENSEMBL_GTF_FILENAME` via
  `urllib.request.urlretrieve`; skips the download if the target already exists unless
  `--overwrite` is passed.
- **Drift probe:** `fetch_fingerprint` in
  `hvantk/skills/ensembl_gene/structure/drift_probe.py`, `PROBE_VERSION = 2`. HEADs
  `ENSEMBL_GTF_URL` (no payload download -- the docstring notes the GTF is ~64 MB; note `catalog/datasets.json` records `size_bytes: 800000000` for the same file, so one of the two is wrong and neither should be trusted for capacity planning) and
  fingerprints `ETag` + `Content-Length` under `headers`, records `Last-Modified` as
  `source_version`, and puts `{release, url}` under `extras`. Raises `DriftProbeError` on
  an HTTP failure or if the response has neither `ETag` nor `Content-Length`.
- **Release pin:** `hvantk/resources/ensembl_release.py`, also consumed by
  `hvantk.algorithms.ptm.constants` and by the separate parser
  `hvantk.algorithms.ptm.mapper.parse_ensembl_gtf` (see § 1).
- **Catalog:** `hvantk/skills/ensembl_gene/catalog/datasets.json`.
- **CLI:**
  `hvantk reprocess ensembl-gene:structure --raw-dir <dir> --output <out>.ht [--plugin-arg protein_coding_only=true]`;
  standalone download via `hvantk download ensembl-structure --raw-dir <dir> [--overwrite]`.
  (The manifest's `cli:` command is `ensembl-structure-download`; the loader strips the
  `-download` suffix and binds it under the `download` group, so the suffixed form is
  not a command.)
- **Tests:** `hvantk/skills/ensembl_gene/structure/tests/` -- `test_builder.py` (Hail
  snapshot round-trip + key-uniqueness + `protein_coding_only` filter + raw-directory
  handling), `test_parse.py` (Hail-free parser unit tests), `test_resolve_path.py`
  (Hail-free `_resolve_gtf_path` unit tests), `test_drift_probe.py` (offline, mocked
  drift-probe tests via `requests_mock`).

## 7. Workflow steps

When invoked to build, refresh, or extend the `ensembl-gene:structure` table:

1. **Confirm or fetch the raw GTF.** `hvantk download ensembl-structure --raw-dir <dir>`
   (or let `hvantk reprocess` invoke `download_dataset` itself). This is a no-op if
   `raw_dir/Homo_sapiens.GRCh38.113.gtf.gz` already exists, unless `--overwrite` is
   passed.
2. **Build the table.**
   `hvantk reprocess ensembl-gene:structure --raw-dir <dir> --output <out>.ht [--plugin-arg protein_coding_only=true]`,
   or call `build_ensembl_gene_structure(parsed_input, ctx, protein_coding_only=...)`
   directly. `parsed_input` may be either the raw directory (as `reprocess` passes it,
   per § 4) or a direct GTF file path, gzipped or plain -- both are handled by
   `_resolve_gtf_path`.
3. **Sanity-check the output.** Confirm the table is keyed by `gene_id` with
   `ht.count() == ht.distinct().count()` (see `test_builder.py::test_gene_id_is_the_key_and_is_unique`),
   and that a known MANE-annotated gene has a non-empty `mane_select`.
4. **Run the dataset's own tests:**
   `pytest hvantk/skills/ensembl_gene/structure/tests -m hail` (§ 9). Note
   `test_parse.py` and `test_resolve_path.py` are Hail-free and also collect under a plain
   `pytest` run of the default suite, since they carry no `@pytest.mark.hail`.
5. **If anything about the release changed, verify pin consistency:**
   `pytest hvantk/skills/ensembl_gene/tests/test_release_pin.py` and
   `pytest hvantk/tests/test_ensembl_release_pin.py` (§ 2, § 8).

## 8. Update playbook

Triggered by Ensembl publishing a new release, or by the pinned GTF file being re-issued
upstream under the same release.

1. **Run the drift probe:** `hvantk drift ensembl-gene:structure`. It HEADs
   `ENSEMBL_GTF_URL` and compares `ETag`/`Content-Length` against
   `structure/tests/drift_fingerprint.json`. By design it never reports "a newer Ensembl
   release exists" as drift -- only movement of the *pinned* file's own validators (see
   `drift_probe.py`'s module docstring: "The release is a repo-side pin, not an upstream
   fact"). Deciding to move to a new release is a separate, manual step below.
2. **To adopt a new release,** bump `ENSEMBL_RELEASE` in
   `hvantk/resources/ensembl_release.py`. `ENSEMBL_GTF_FILENAME` and `ENSEMBL_GTF_URL`
   are both f-strings derived from it, so this one edit propagates to every consumer.
3. **Update the catalog entry** in
   `hvantk/skills/ensembl_gene/catalog/datasets.json` -- `accession` must become
   `Ensembl_v<release>` and the `.gtf.gz` file entry's `path` must match the new
   `ENSEMBL_GTF_FILENAME`; `test_release_pin.py::test_catalog_declares_the_same_release`
   enforces this.
4. **Re-download and rebuild:**
   `hvantk download ensembl-structure --raw-dir <dir> --overwrite`, then
   `hvantk reprocess ensembl-gene:structure --raw-dir <dir> --output <out>.ht` end to end.
5. **Regenerate the committed drift baseline,**
   `structure/tests/drift_fingerprint.json`, once the change is validated -- not silently
   in the same PR as a behavioural change (`_conventions` § 12).
6. **Re-run both pin tests:**
   `hvantk/skills/ensembl_gene/tests/test_release_pin.py` and
   `hvantk/tests/test_ensembl_release_pin.py`. The second fails if
   `hvantk.algorithms.ptm.constants` was not updated to match -- that cross-layer
   assertion is the intended guard against the two pins drifting apart.
7. **If the fixture needs new coverage** (e.g. a new MANE/tie-break shape), extend
   `structure/tests/testdata/raw/ensembl-structure/mini.gtf` and regenerate snapshots:
   `pytest hvantk/skills/ensembl_gene/structure/tests/test_builder.py -m hail --regenerate-snapshots`.

## 9. Validation contract

Declared in `plugin.yaml`'s `datasets[0].tests` block, paths relative to the provider
folder (`hvantk/skills/ensembl_gene/`):

- `fixture`: `structure/tests/testdata/raw/ensembl-structure` -- i.e.
  `hvantk/skills/ensembl_gene/structure/tests/testdata/raw/ensembl-structure/mini.gtf`.
  **This is a plugin-local (dataset-local) fixture, not the shared repo-level form.**
  `_conventions` § 9 currently lists `ensembl_gene` among providers using the shared
  `../../tests/testdata/raw/<dataset>` fixture path; verified against the manifest and
  the test files above, that is stale for the `structure` dataset -- `plugin.yaml`
  declares, and `test_builder.py`/`test_parse.py` read, a fixture committed inside
  `structure/tests/` itself.
- `schema_snapshot`: `structure/tests/snapshots/schema.json` -- seeded.
- `row_snapshot`: `structure/tests/snapshots/sample_rows.json` -- seeded, 3 rows keyed by
  `gene_id` (one per fixture gene).
- `drift_fingerprint`: `structure/tests/drift_fingerprint.json` -- seeded
  (`probe_version: 2`, `extras.release: "113"`).
- `command`: `pytest hvantk/skills/ensembl_gene/structure/tests -m hail`.

Not part of `plugin.yaml`'s `tests:` block, but load-bearing companions to this
dataset's contract (§ 2, § 8): `hvantk/skills/ensembl_gene/tests/test_release_pin.py`
(plugin-local, Hail-free) and `hvantk/tests/test_ensembl_release_pin.py` (repo-level,
Hail-free, cross-checks `hvantk.algorithms.ptm.constants`).

> **Snapshot status:** seeded. `structure/tests/snapshots/schema.json` and
> `structure/tests/snapshots/sample_rows.json` are committed and
> `test_structure_snapshot_round_trip` asserts the build against them. Regenerate after
> an intentional schema change with
> `pytest hvantk/skills/ensembl_gene/structure/tests/test_builder.py -m hail --regenerate-snapshots`,
> then commit the result.
