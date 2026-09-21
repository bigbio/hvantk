---
name: hvantk:resource-ensembl-gene
description: Ensembl gene-model provider -- ships the ensembl-gene:structure dataset, the canonical per-gene structural table parsed from the pinned-release GTF.
status: provisional
backend: hail
domain: mapping
---

# ensembl-gene

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

This is the provider-level file. It exists so every `SKILL.md` under `hvantk/skills/`
carries the same nine-section contract, but each section here is brief and points at the
one dataset's file for the substance -- see `structure/SKILL.md`, which is authoritative.

## 1. Status & scope

Provisional. This provider ships exactly one dataset today, declared in
`hvantk/skills/ensembl_gene/plugin.yaml`: `ensembl-gene:structure` -- the canonical
per-gene Ensembl table (gene ID, name, biotype, coordinates, CDS length, coding-exon
count, transcript count, MANE Select) parsed from the pinned-release GTF, and the gene
spine the annotation pipeline builds on. Full status/scope/out-of-scope notes are in
`structure/SKILL.md` § 1.

## 2. Source identity

Provider: Ensembl (<https://www.ensembl.org/info/data/ftp/index.html>). Catalog entry:
`hvantk/skills/ensembl_gene/catalog/datasets.json` (`accession: Ensembl_v113`,
`data_source: "Ensembl"`). The release pin (`hvantk/resources/ensembl_release.py`) is
shared by this provider's dataset and by `hvantk.algorithms.ptm`; the full pin-identity
and machine-checked consistency contract is in `structure/SKILL.md` § 2.

## 3. Backend choice + reasoning

`backend: hail`, `domain: mapping`, declared per-dataset in `plugin.yaml` (there is one
dataset, so provider-level and dataset-level agree). Reasoning for choosing Hail is
dataset-specific and lives in `structure/SKILL.md` § 3.

## 4. Raw format & gotchas

GTF parsing gotchas (feature-row filtering, regex-based attribute extraction,
gene/transcript version stripping, MANE-vs-longest-CDS tie-break, and the
directory-vs-file input each `reprocess` build hands the builder) are all
dataset-specific. See `structure/SKILL.md` § 4 -- not restated here.

## 5. Output contract

Every dataset under this provider returns an `AnnotationTable` (Hail Table). The one
current dataset's exact schema, key, and row contract are in `structure/SKILL.md` § 5.

## 6. hvantk integration points

- **Manifest:** `hvantk/skills/ensembl_gene/plugin.yaml` (`api_version: 2`, provider
  name `ensembl-gene`).
- **Catalog:** `hvantk/skills/ensembl_gene/catalog/datasets.json`.
- **Dataset folder:** `hvantk/skills/ensembl_gene/structure/` -- builder, parser,
  downloader CLI, and drift probe; see `structure/SKILL.md` § 6 for each callable.
- **Provider-wide pin-consistency tests** (apply across anything under this provider
  that reads the release pin): `hvantk/skills/ensembl_gene/tests/test_release_pin.py`
  (plugin-local) and `hvantk/tests/test_ensembl_release_pin.py` (repo-level; cross-checks
  `hvantk.algorithms.ptm.constants`).

## 7. Workflow steps

`hvantk reprocess ensembl-gene:structure --raw-dir <dir> --output <out>.ht [--plugin-arg protein_coding_only=true]`
builds the sole dataset end to end. The full step-by-step (download, build,
sanity-check, test) is in `structure/SKILL.md` § 7; this provider adds no cross-dataset
orchestration on top since it ships only one dataset.

## 8. Update playbook

Triggered by an Ensembl release bump or an upstream re-issue of the pinned GTF. The
playbook (bump `ENSEMBL_RELEASE`, update the catalog accession, re-download, re-run the
pin tests, regenerate the drift fingerprint) is entirely dataset-scoped today -- see
`structure/SKILL.md` § 8 for the full sequence.

## 9. Validation contract

Declared in `plugin.yaml`'s `datasets[0].tests` block; the exact `fixture` /
`schema_snapshot` / `row_snapshot` / `drift_fingerprint` / `command` paths, and the
note that the fixture is plugin-local (not the shared repo-level form
`_conventions` § 9 otherwise implies for this provider), are in `structure/SKILL.md` § 9.
No provider-level test target exists beyond that dataset's own
`pytest hvantk/skills/ensembl_gene/structure/tests -m hail` plus the two pin-consistency
tests named in § 6.
