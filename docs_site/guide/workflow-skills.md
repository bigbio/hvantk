# Workflow skills

hvantk ships **workflow skills**: instructions a coding agent follows to run the toolkit
for you, on a laptop or on a Slurm cluster, with the checks a careful operator runs.
The agent orchestrates; every command it executes is a deterministic hvantk command,
and no language model sits in the data path.

They live in the repository's top-level `skills/` folder, one folder per skill, in the
Agent Skills layout. The catalog and the install steps are in
[`skills/README.md`](https://github.com/bigbio/hvantk/tree/main/skills); the contract
they follow is `skills/workflow-conventions/SKILL.md`.

## What a skill gives you

- A **procedure**: the exact commands, in order, with the template scripts to run them.
- **Preflight** that refuses to guess: the checkout commit, the image digest and every
  cluster-specific value are confirmed first, from a site profile you fill in once.
- **Traps**: the failure modes we met, each with its cause and fix.
- **Verification**: what to check on a real run and a receipt tied to the output's
  provenance sidecar.
- A **measured task**: each skill ships an `eval/` folder with a public-data task, a
  deterministic check and the reference values its positive control produced, so the
  skill can be tested, by us and by you.

## Getting started

Fill in the site profile on the cluster, in the `$HVANTK_REPO` checkout:

```bash
cp skills/_site/site-profile.example.sh skills/_site/site-profile.local.sh   # gitignored
# fill it in: see skills/_site/SITE-PROFILE.md
export HVANTK_SITE_PROFILE=/path/to/hvantk/skills/_site/site-profile.local.sh
bash skills/_site/check-site-profile.sh
```

No image yet on this site? `check-site-profile.sh` otherwise refuses because
`HVANTK_SIF` points nowhere; run the `hvantk:hpc-container` skill first, or pass it
`--allow-missing-image` while bootstrapping.

Then install the skills into your agent (Claude Code: `/plugin install hvantk --marketplace
bigbio/hvantk`, or on older releases `/plugin marketplace add bigbio/hvantk` then
`/plugin install hvantk@hvantk`; other harnesses: point them at `skills/`) and ask
it, for example, to build a dataset on the cluster. The shipped skills cover running any
Hail-backed hvantk command as a single-node Slurm job from the container image, keeping a
cluster clone in step with GitHub, building and promoting the image, joint genotyping,
delivering a sample subset of a callset, ancestry inference, building many datasets as
one array job with their outputs verified, triaging a drift signal (repairing the
plugin, rebuilding, and recording the fingerprint and the rebuild ledger), building the
gene index, running a cohort's gene burden and attaching it, ranking a cohort's genes and
evaluating the ranking against a label set, and the molecular QTL cascade with
colocalisation; the catalog in `skills/README.md` is the current list. See the [HPC migration guide](hpc-migration.md)
for the cluster setup those skills assume.

## The two kinds of skill

The per-provider `SKILL.md` files under `hvantk/skills/<provider>/` are **resource
specs**: they tell an agent how to author a data-source plugin against the plugin
contract. **Workflow skills** tell it how to run the toolkit. They reference each other
and do not overlap.
