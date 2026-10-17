# hvantk workflow skills

Skills that let a coding agent drive an hvantk analysis for you: on a laptop or on a
Slurm cluster, from the container image, with the checks a careful operator would run.
The agent orchestrates; every command it runs is a deterministic hvantk command, and no
language model sits in the data path.

There are two kinds of skill in this repository. **Resource specs** under
`hvantk/skills/<provider>/SKILL.md` tell an agent how to author a data-source plugin.
**Workflow skills**, here, tell it how to run the toolkit. Read
`workflow-conventions/SKILL.md` for the contract they follow.

## Catalog

| skill | use it when |
| --- | --- |
| `hvantk:hpc-hail-job` | you want to run any Hail-backed hvantk command (`reprocess`, `hgc`, `annotate`, `cohort`, `qtlcascade`) on a Slurm cluster from the container image: sizing, scratch, submission, waiting, triage, verification |
| `hvantk:hpc-sync` | a cluster clone or bare remote must be brought in step with GitHub without losing cluster-side commits, branches must be pushed to the cluster remote by explicit refspec, or a job ran a commit you did not review |
| `hvantk:hpc-container` | the container image must be built from the lock file without root, validated with the install check, recorded with provenance, or promoted as the default; or a job runs old dependencies although the code is current |
| `hvantk:joint-genotyping` | single-sample gVCFs must become a joint callset (VDS, MatrixTable, QC, cohort VCF) with the HGC pipeline, for a few samples in one job or thousands per contig |
| `hvantk:callset-subset-export` | a sample subset of a joint callset must be delivered as a verified VCF with checksums and a manifest, by subsetting the lossless VDS, per contig for large callsets |
| `hvantk:ancestry-inference` | samples must be assigned a continental ancestry against a labelled reference built from 1000 Genomes, with the PCA and random-forest pipeline as a job |
| `hvantk:dataset-build-at-scale` | one or many `hvantk reprocess` datasets must be built from staged raw files as Slurm jobs — fanned out as one array job, each output verified against its plugin's committed schema contract |
| `hvantk:drift-triage` | `hvantk drift` or a drift-bot pull request reports a dataset drifted or its probe failed: classify the signal, repair the plugin or accept the fingerprint, rebuild and verify against the committed snapshots, and record the fingerprint, the rebuild-ledger row and the rebuild |
| `hvantk:gene-index` | the gene index that gene-level annotations, cohort tables and rankings join onto must be built from the Ensembl gene-structure and HGNC datasets (`hvantk annotate spine`), its gene count and HGNC mapping rate checked, and annotation axes added onto it |
| `hvantk:cohort-attach-burden` | a case/control cohort must become a gene-level table: the gene burden test on a cleaned MatrixTable, the cohort manifest, validation, and attachment to the gene index or the annotation matrix |
| `hvantk:rerank-evaluate` | a cohort's genes must be ranked by multi-omic evidence with `hvantk rerank` and evaluated against a label set you choose (ClinGen or GenCC), with the identifiers aligned, paralogue blocks and the case/control audit |
| `hvantk:qtl-cascade` | eQTL and pQTL tables must be built and run through the molecular QTL cascade with colocalisation and its report for one or more tissues, or a GWAS locus colocalised with eQTLs |

## Install

**Claude Code**, from the repository:

```text
/plugin install hvantk --marketplace bigbio/hvantk
```

or, on older Claude Code releases, in two steps: `/plugin marketplace add bigbio/hvantk`
then `/plugin install hvantk@hvantk`.

or from a checkout, for development: `claude --plugin-dir /path/to/hvantk`.

**Other agent harnesses** that read the Agent Skills layout: point them at this
`skills/` folder, or symlink it into their skills path. hvantk's own keys sit under
`metadata`, as the layout requires; the two Claude Code keys (`user-invocable`,
`argument-hint`) are ignored by other harnesses.

## Before the first run: the site profile

Every cluster-specific value lives in one file you fill in once, on the cluster, in
the `$HVANTK_REPO` checkout:

```bash
cp skills/_site/site-profile.example.sh skills/_site/site-profile.local.sh   # gitignored
# fill it in (see skills/_site/SITE-PROFILE.md), then:
export HVANTK_SITE_PROFILE=/path/to/hvantk/skills/_site/site-profile.local.sh
bash skills/_site/check-site-profile.sh
```

A skill that finds a value missing stops and asks. It never guesses a partition, a
path or an image.

No image yet on this site? `check-site-profile.sh` otherwise refuses because
`HVANTK_SIF` points nowhere; run the `hvantk:hpc-container` skill first, or pass it
`--allow-missing-image` while bootstrapping.

## How a skill is measured

Each skill ships an `eval/` folder: a task an agent can be handed cold, a script that
stages public inputs, a positive control that runs the procedure with no agent, a
deterministic check, and the reference values that control produced. The contract is in
`workflow-conventions/SKILL.md` section 5. The evaluation harness that runs agents against these
tasks is separate from this repository; the tasks and checks are public so anyone with a
cluster can repeat them.
