---
name: hvantk:resource-clinvar
description: Onboard, build, or update the ClinVar resource for hvantk
status: provisional
backend: hail
domain: variants
---

# ClinVar resource skill

## 1. Status & scope

This skill covers BUILD and UPDATE of the ClinVar Hail Table. It does NOT cover download (see `hvantk/skills/clinvar/cli.py`) or downstream analysis.

## 2. Source identity

Provider metadata (URL, version cadence, license, citation) lives in the ClinVar entry inside `hvantk/resources/registry/genomics/datasets.json` (find by `accession: "ClinVar_latest"`), surfaced via `hvantk catalog show ClinVar_latest`. Read that file; do not restate.

Stable note (not in catalog): ClinVar releases monthly. New INFO fields are uncommon but do occur (e.g., the addition of `ONCOGENICITY` in the oncogenicity supplement).

## 3. Backend choice + reasoning

Hail Table keyed by `(locus, alleles)`. Reasoning: ClinVar is variant-keyed and joins with other variant tables (gnomAD, dbNSFP) downstream. A single Hail Table is the natural fit; nothing about ClinVar requires anndata or pandas.

## 4. Raw format & gotchas

- File format: bgzipped VCF (`.vcf.bgz` with `.tbi` index).
- Reference genome: GRCh38 by default. The `contig_recoding()` helper from `hvantk/tables/table_builders.py` maps numeric contigs to `chr*` names.
- Import: use `hl.import_vcf(path, force=True, reference_genome=..., contig_recoding=..., skip_invalid_loci=True)`. The `force=True` flag is required because ClinVar VCFs contain non-standard headers.
- Keying: after `.rows()`, repartition (typical: 100) and key by `(locus, alleles)`.
- INFO fields commonly used: `CLNSIG`, `CLNREVSTAT`, `CLNDN`, `CLNDISDB`, `RS`, `MC`, `GENEINFO`. The full list comes from the VCF header — read it, do not assume.
- Encoding gotchas:
  - Spaces in INFO values are encoded as `_` (e.g., `Likely_pathogenic`).
  - Multi-value INFO fields use `,` or `|` depending on the field.
  - `CLNDISDB` uses `,` between databases and `|` between IDs within a database.
- TSV export gotcha: when `export_tsv=True`, ClinVar's TSV is **flattened** before export so nested struct fields become flat columns. Do NOT pass `export_tsv` straight through to `_create_table_base` — instead, force `export_tsv=False` in the helper call and run `<table>.flatten().export(f"{output_path}.tsv.bgz")` after the helper returns. Other builders may forward `export_tsv` directly; ClinVar is special because of the deeply nested INFO struct.

## 5. Output contract

Hail Table at `<output_path>.ht`. Schema is defined by `hvantk/tests/snapshots/clinvar/schema.json` (canonical). Human summary: keyed by `(locus, alleles)`; row contains `rsid`, `qual`, `filters`, and an `info` struct with the ClinVar-specific INFO fields parsed by `hl.import_vcf`.

## 6. hvantk integration points

- Builder: `create_clinvar_tb` in `hvantk/skills/clinvar/builder.py`
- Downloader CLI: `clinvar_downloader` in `hvantk/skills/clinvar/cli.py` (registered as `hvantk download clinvar`)
- Dataset class: `ClinVarDataset` in `hvantk/skills/clinvar/shared/datasets.py`
- Build CLI: `hvantk mktable clinvar` in `hvantk/tools/build/make_table_cli.py`
- Plugin manifest: `hvantk/skills/clinvar/plugin.yaml` (drives loader registration; compound dataset key `clinvar:variants`)
- Test: `hvantk/skills/clinvar/tests/test_builder.py`

Read the existing files at these paths as ground truth for shape. This skill does not restate code.

## 7. Workflow steps

When invoked to build or update:

1. Verify Hail is available (defer to the SessionStart hook).
2. Confirm input is a ClinVar VCF: read the `#CHROM` header line of the input file.
3. Use `_create_table_base()` from `hvantk/tables/table_builders.py` with:
   - `import_func = lambda: hl.import_vcf(input_path, force=True, reference_genome=reference_genome, contig_recoding=contig_recoding(), skip_invalid_loci=True).rows()`
   - `transform_func = lambda ht: ht.repartition(100).key_by("locus", "alleles")`
4. Apply the parsing rules from § 4 (force, contig_recoding, skip_invalid_loci).
5. Preserve the flatten-before-export branch from § 4 — pass `export_tsv=False` to `_create_table_base` and run `<table>.flatten().export(...)` after. The round-trip test does NOT exercise `export_tsv=True`, so a regression here would land silently.
6. Run validation: `pytest hvantk/skills/clinvar/tests/test_builder.py -m hail`.
7. Report: schema diff, sample-row diff, test pass/fail.

## 8. Update playbook

When ClinVar releases a new monthly version:

1. Fetch the new release: `hvantk download clinvar --output-dir /tmp/clinvar`.
2. Regenerate the fixture: extract a small representative slice (mix of CLNSIG values, multi-allelic site, multi-CLNDN row). The current pilot fixture is a single-chromosome slice (`hvantk/skills/clinvar/tests/testdata/raw/clinvar/clinvar_20220403_chr20.vcf.bgz`) — this is adequate because ClinVar parsing is INFO-field driven, not chromosome-dependent. Replace the file (re-bgzip if needed). A `.tbi` index is optional; `hl.import_vcf(force=True)` reads `.bgz` directly.
3. Run snapshot regeneration: `pytest hvantk/skills/clinvar/tests/test_builder.py -m hail --regenerate-snapshots`.
4. Inspect the snapshot diff:
   - **Expected diff** (new INFO field, additional CLNSIG value): commit the regenerated snapshots with explanation.
   - **Unexpected diff** (schema regression, missing field): STOP. Investigate before committing.
5. If a parsing gotcha was introduced (e.g., new encoding), add it to § 4.
6. Open PR; reviewer checks the snapshot diff narrative.

## 9. Validation contract

- `fixture`: `hvantk/skills/clinvar/tests/testdata/raw/clinvar/clinvar_20220403_chr20.vcf.bgz`
- `schema_snapshot`: `hvantk/skills/clinvar/tests/snapshots/schema.json`
- `row_snapshot`: `hvantk/skills/clinvar/tests/snapshots/sample_rows.json`
- `test_command`: `pytest hvantk/skills/clinvar/tests/test_builder.py -m hail`
