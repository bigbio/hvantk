---
name: ml-features
description: 'Use when building a labelled machine-learning training set of ClinVar variants with hvantk training-set build on a Slurm cluster: pathogenic/benign labels with review stars, dbNSFP scores and allele frequencies, protein interface and PTM flags, gene-level features from the annotation matrix, gene-list flags, and a leakage tier (L1-L4) for every feature column. Also when choosing training columns free of label or allele-frequency leakage (circularity, trained_on), when a feature column comes out as text, or when the build stops with "--input does not match the spec", "gene matrix lacks column(s)", "is not in the gene matrix''s manifest", "that the matrix reducer did not produce" or "is not a field of".'
user-invocable: true
argument-hint: "<spec.yaml> [--gene-list COLUMN=FILE] [--out PREFIX]"
metadata:
  stage: analysis
  status: provisional
  scale: node
---

# Build a labelled, leakage-tiered ML training set

Read the `hvantk:workflow-conventions` skill (`skills/workflow-conventions/SKILL.md`)
first; this skill assumes its site-profile rule, its section contract and its
measurement contract. Jobs run through `hpc-hail-job`.

## 1. Scope and outcome

A training set built by `hvantk training-set build` from one feature-spec YAML, as four
files under one output prefix:

- `<prefix>.ht`, keyed `locus`, `alleles`, with a `.provenance.json` sidecar;
- `<prefix>.tsv.gz`, the same rows as a deterministic gzipped TSV (byte-identical when
  rebuilt from the same inputs);
- `<prefix>.dictionary.tsv`, one row per column: whether it is a training feature, its
  tier, level, source and source version, non-null and true rates;
- `<prefix>.manifest.json`: the row funnel, label balance (overall and per gene list),
  review stars, how each row's gene was chosen, tier counts, excluded columns, gene-list
  mapping, and every input's provenance;

plus a receipt from `templates/describe-training-set.py` that checks the four agree and
binds them to the `.ht` sidecar's and the TSV's SHA-256.

Rows are ClinVar GRCh38 SNVs with the spec's consequence (missense by default), labelled
1 or 0 from the primary CLNSIG, with at least the spec's review stars. The shipped spec,
`examples/training_sets/clinvar_heart.yaml`, is the worked example: 44 dbNSFP rankscores,
three gnomAD allele frequencies, INSIDER interface and UniProt PTM flags, six gene axes
(constraint, GeVIR, interaction partners, three heart expression atlases) and a CHD
gene-list flag.

Use it when a model needs ClinVar-labelled variants with predictor scores and gene
features, or when someone asks which columns can be trained on without leakage. Not this
skill: building the gene index (`gene-index`); ranking genes (`rerank-evaluate`); labels
other than ClinVar's (unsupported by the command).

## 2. Preflight

1. Site profile:
   `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`;
   stop and ask for any missing value. Every later command that reads a profile value
   starts with the same `source "$HVANTK_SITE_PROFILE" &&`.
2. Checkout and its commit: `git -C "$HVANTK_REPO" rev-parse HEAD` and
   `git -C "$HVANTK_REPO" status --porcelain`; compare with the commit the user expects.
   The checkout must have the command:
   `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- hvantk training-set build --help`.
3. Image and its digest: `readlink -f "$HVANTK_SIF"` and
   `sha256sum "$(readlink -f "$HVANTK_SIF")"`.
4. Inputs and outputs: every built input, the gene index and the output prefix are
   readable or writable from compute nodes, under `$HVANTK_WORK` (or bound through
   `HVANTK_BIND` in the profile file).
5. The spec: the shipped one or an edited copy (`references/spec-and-tiers.md`). Every
   `--input` name it needs:
   `grep -nE '^\s+(- )?(input|hgnc_input):' <spec.yaml>`. Every `gene_lists` column needs a
   `--gene-list COLUMN=FILE` (symbols separated by whitespace); ask the user for each
   file, never guess one.
6. The ClinVar release: a dated release (`YYYYMMDD`), decided now and written into the
   run notes. `latest` moves every week, and a different release is a different dataset.
7. Inputs that already exist: for each, the sidecar's `builder_commit` against
   `git -C "$HVANTK_REPO" rev-parse HEAD`, and for dbNSFP the type of its
   `*_rankscore` fields, which must be `float64` (`references/inputs.md` has the
   check). A dbNSFP table with `str` rankscores is rebuilt, never patched with `cast:`.
   A `builder_commit` of `hvantk <version>` means the job's node had no `git` (its log
   says `commit: unknown`): the commit is then only in that build's run notes.
8. The gene index: an `annotate spine` output from the `gene-index` skill, built with the
   same `hgnc:lookup` table the spec's `hgnc` input will name.

## 3. Procedure

`RAW` and `OUT` below are absolute directories under `$HVANTK_WORK`; `<run>` is the run id.
Copy the spec to `$OUT/spec.yaml` and each gene list next to it first; the commands read
them there.

1. **Stage the raw files on the login node**, through the image, one plugin downloader
   each, ClinVar at the pinned release (`--verify-md5` checks it against NCBI's
   checksum; keep the VCF, the only lasting record of which release was read):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh" -- \
     hvantk download clinvar --version YYYYMMDD --verify-md5 --output-dir "$RAW/clinvar"
   ```
   `references/inputs.md` lists every input of the shipped spec with its downloader, the
   two that have none (dbNSFP, the INSIDER variants BED) and how to get them.
2. **Build every input dataset as a job**, the same shape for each:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-clinvar --cpus 8 --mem-gb 32 --time 0-01:00:00 --wait -- \
     hvantk reprocess clinvar:variants --raw-dir "$RAW/clinvar" --output "$OUT/inputs/clinvar.ht" \
       --skip-download --no-check-drift
   ```
   or fan the small ones out as one array with `dataset-build-at-scale`. dbNSFP is far
   bigger than the rest: its own whole-node job
   (`skills/dataset-build-at-scale/references/sizing-by-dataset.md`), built once and reused.
3. **Summarise expression for every `layer1` axis with a `matrix:` block** (the shipped
   spec has three); the commands are in `references/inputs.md`. Each summary must be
   grouped on the field the spec's column names come from.
4. **Build the gene matrix from the same spec**: one `annotate prepare` job per `layer1`
   axis, then one `annotate compose` job with exactly one `--prepared AXIS=PATH` per axis
   (`skills/gene-index/references/annotation-axes.md` has both commands and the `--hgnc`
   rule):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-prep-gevir --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
     hvantk annotate prepare --spec "$OUT/spec.yaml" --axis gevir --input "$OUT/inputs/gevir.ht" \
       --spine "$OUT/gene-index.ht" --output "$OUT/prepared/gevir.ht"
   ```
   `compose` writes `<output>.manifest.json` next to the matrix; the training-set build
   reads it.
5. **Build the training set as a job**:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-training-set --cpus 8 --mem-gb 48 --time 0-01:00:00 --wait -- \
     hvantk training-set build --spec "$OUT/spec.yaml" \
       --input clinvar="$OUT/inputs/clinvar.ht" --input dbnsfp="$OUT/inputs/dbnsfp.ht" \
       --input insider_variants="$OUT/inputs/insider_variants.ht" --input ptm="$OUT/inputs/ptm_sites.ht" \
       --input hgnc="$OUT/inputs/hgnc-lookup.ht" --gene-matrix "$OUT/gene_matrix.ht" \
       --gene-list is_chd_gene="$OUT/chd_genes.txt" --output "$OUT/training_set/clinvar_heart"
   ```
   One `--input` per name from Preflight step 5, no more and no fewer. Spec, input and
   gene-matrix mistakes stop the command before any Hail job runs (Traps).
6. **Describe it and write the receipt**, as a job:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" \
     --name <run>-describe --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python "$HVANTK_REPO/skills/ml-features/templates/describe-training-set.py" \
       --prefix "$OUT/training_set/clinvar_heart" --receipt "$OUT/training_set/clinvar_heart.receipt.json"
   ```
   It exits 1 when the four outputs disagree or a feature column is typed `str`.
7. **Run notes**: `commit:`, `container:`, one `job:` line per job, the manifest's
   `spec.sha256`, and the ClinVar file itself: its name, its `sha256sum` and its
   `##fileDate` line (`zgrep -m1 '^##fileDate=' RAW/clinvar/<file>.vcf.gz`). The
   dataset does not record which release it read: its `source_fingerprint` is the
   upstream's state when the build ran (or `probe-unavailable` offline), not the file.

## 4. Scale-up

The shipped spec, on ClinVar of November 2026 and dbNSFP 4.9a: 4,617,349 ClinVar records,
2,595,328 missense SNVs, 227,513 labelled, 210,867 with at least one star; 312 columns,
301 of them features. With the inputs and the gene matrix already built, one
`training-set build` job took 5 minutes, with a peak of 21 GB resident memory on a
32-core allocation that stayed almost idle (2% CPU efficiency): the request in step 5
leaves headroom; read `seff` after the first run and adjust. The TSV is written through
an uncompressed export copied to node-local disk (850 MB for that build, 59 MB
compressed): budget a few GB of node-local scratch.

The cost is in the inputs, not in the build: dbNSFP's academic release is about 50 GB to
download and a table of about 70 GB once built, built once on a whole node and reused;
the expression summaries read GB-sized matrices; the rest take minutes. Rebuild the training set when
the ClinVar release, any input, the gene matrix or the spec changes; rebuild the gene
matrix when the spec's `layer1` or the gene index changes.

## 5. Traps

- Feature columns typed `str`, values like `"0.71234"` in the TSV. Cause: a dbNSFP table
  built before the builder parsed `*_rankscore` as numbers; `training-set build` passes
  a string column through without an error. Fix: rebuild `dbnsfp:variants` with this
  checkout; the receipt's `string_features` must be empty.
- `--input does not match the spec (missing for: ...; not declared in the spec: ...)`
  (or the same for `--gene-list`). Cause: the names on the command line differ from the
  spec's `input:`, `hgnc_input:` and `gene_lists` keys. Fix: pass exactly those names.
- `axis '<axis>' is not in the gene matrix's manifest; was the matrix composed from this
  spec?`. Cause: the matrix was composed from another spec, or before an axis was added.
  Fix: prepare and compose again from this spec. The build checks only that each axis
  is there: a matrix composed before an axis's `columns`, `collapse` or `stats` changed
  passes with the old values, so prepare and compose again after any `layer1` edit.
- `gene matrix lacks column(s): ...`. Cause: the matrix and its compose manifest
  disagree (a `--gene-matrix-manifest` from another compose), or the gene index lacks a
  column the spec's `gene.spine_columns` names. Fix: pass the manifest written next to
  the matrix; build the gene index with `annotate spine`.
- `entry '<source>' declares column(s) [...] that the matrix reducer did not produce;
  it emitted [...]`, from `annotate prepare`. Cause: an expression summary grouped on
  another field, or an atlas whose group labels differ from the spec's. Fix: summarize
  again with the `--group-by` fields `references/inputs.md` gives for that axis.
- `gene-matrix manifest '<path>' not found`. Cause: the compose manifest is not next to
  the matrix. Fix: pass `--gene-matrix-manifest`.
- `field '<path>': '<segment>' is not a field of ...`. Cause: a spec column, a
  `gene_field`, or a codon source's `chrom`/`codon_start`/`codon_end` names a field the
  input does not have. Fix: correct the spec against the input's schema.
- A column missing from the dataset and listed under `excluded`. Cause: neither a
  `trained_on` declaration nor a `tier:` (every `layer1` axis needs one; an unregistered
  source logs `its columns need a tier in the spec`). Fix: declare a tier, or leave the
  column out on purpose.
- A `tier: L1` declaration that comes out L2 or L3. Cause: by design, a declaration only
  raises a tier; the plugin's `trained_on` (ClinVar, HGMD, OMIM, population frequency)
  sets the floor. Fix: none; that column is not clean.
- Two runs of the same spec giving different datasets. Cause: ClinVar `latest` (or
  another moving download) changed between them. Fix: pin the release; compare the run
  notes' ClinVar file, `sha256sum` and `##fileDate`, and the receipts'
  `funnel.clinvar_rows`. The manifests' `input_provenance` does not tell releases apart.
- Variants at the wrong positions, joined to the wrong dbNSFP rows, with no error.
  Cause: a GRCh37 ClinVar VCF (`hvantk download clinvar --genome-build GRCh37`); the
  builder reads every VCF as GRCh38 and skips the loci that do not fit. Fix: GRCh38
  only, for every input keyed by position (gene-keyed axes such as gnomAD constraint
  v2.1.1 join on `gene_id` and are unaffected).
- `clinvar_clnsig` and `clinvar_review_status` with commas inside one value. Cause:
  `import_vcf` splits these INFO fields at commas and the build joins them back, so the
  value is the one ClinVar released. Fix: none; do not split them again downstream.
- `gene_id` and `gene_symbol` empty on a few rows (`no_gene_id` in the manifest), or a
  multi-gene variant assigned the second listed gene. Cause: none of the variant's
  ClinVar genes maps to an Ensembl id, or the preferring source (dbNSFP's
  `Ensembl_geneid`) assigns that gene. Fix: expected; `gene_rule` counts each case.
- A `--gene-list` symbol that does not map. Cause: an alias, a withdrawn symbol, or a
  gene the gene index does not hold (a non-coding gene on a protein-coding index). Fix:
  it is reported under `gene_lists` in the manifest, not raised; correct the list if
  the count matters.
- `hvantk: command not found`, or a write landing in the wrong checkout. Cause and fix:
  see `workflow-conventions`'s Runtime rules (section 4).

## 6. Verification

Run step 6 on every build and read the receipt (`templates/describe-training-set.py`
prints the same values):

- `consistent: true` and `problems: []`: the `.ht`, the TSV and the manifest have the
  same rows; the TSV header and the dictionary list the `.ht` columns in order; every
  training feature has a tier; the manifest's tier counts equal the dictionary's.
- `string_features: []`.
- `funnel`: each count at most the one before it; `rows` equal to the last.
- `label_counts`: both `0` and `1` present; their balance, overall and per gene list
  (`label_by_gene_list`), is what the user expects for the gene set.
- `tier_counts` and `excluded`: for the shipped spec at dbNSFP 4.9a, L1 265, L2 24,
  L3 8, L4 4, nothing excluded. L1 also depends on the expression summaries' groups
  (each group is a column), so another atlas build changes it. Anything in `excluded`
  must be a column the user chose to leave out.
- `gene_lists`: symbols in, mapped, and genes; `no_gene_id` small (16 of 210,867 rows
  for the shipped spec).
- `inputs`: each input's `builder_commit`, matching Preflight. Which ClinVar release
  was read is in the run notes (step 7), not here.
- `tsv_sha256` and `sidecar_sha256`: the identity of this dataset; a rebuild from the
  same inputs reproduces `tsv_sha256` exactly.

Which columns to train on is read from the dictionary's `training_feature` and `tier`
columns, never from column names: `references/spec-and-tiers.md` has the four training
arms and the selection code.
