# Workflow skills

hvantk ships **workflow skills**: instructions a coding agent follows to run the toolkit
for you, on a laptop or on a Slurm cluster, with the checks a careful operator runs.
The agent orchestrates; every command it executes is a deterministic hvantk command,
and no language model sits in the data path.

They live in the repository's top-level `skills/` folder, one folder per skill, in the
Agent Skills layout. The catalog and the install steps are in
[`skills/README.md`](https://github.com/bigbio/hvantk/tree/main/skills); the contract
they follow is `skills/_conventions/SKILL.md`.

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

```bash
cp skills/_site/site-profile.example.sh skills/_site/site-profile.local.sh   # gitignored
# fill it in: see skills/_site/SITE-PROFILE.md
export HVANTK_SITE_PROFILE=/path/to/hvantk/skills/_site/site-profile.local.sh
bash skills/_site/check-site-profile.sh
```

Then install the skills into your agent (Claude Code: `/plugin marketplace add bigbio/hvantk`
and `/plugin install hvantk@hvantk`; other harnesses: point them at `skills/`) and ask
it, for example, to build a dataset on the cluster. The first shipped skill,
`hpc-hail-job`, runs any Hail-backed hvantk command as a single-node Slurm job from the
container image; see the [HPC migration guide](hpc-migration.md) for how that image is
built.

## The two kinds of skill

The per-provider `SKILL.md` files under `hvantk/skills/<provider>/` are **resource
specs**: they tell an agent how to author a data-source plugin against the plugin
contract. **Workflow skills** tell it how to run the toolkit. They reference each other
and do not overlap.
