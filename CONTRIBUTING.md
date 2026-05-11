# Contributing to hvantk

See the full contributing guide at [docs_site/contributing.md](docs_site/contributing.md).

## Skill maintenance (resource-centric skills)

`hvantk/skills/` contains agent-readable methodology for maintaining stable resources. Pilot scope: ClinVar and UCSC Cell Browser. See `local/planning/2026-05-07-hvantk-resource-skills-design.md` for the full design.

### When you change builder code

If a PR changes a function or path named in any `SKILL.md` — including the `_conventions/` skill — update the affected skill in the same PR. Reviewer will reject PRs where skill drift is obvious.

### Skill status

Each `SKILL.md` carries a `status` frontmatter field:

- `provisional` — work in progress, not yet validated against the round-trip + update tests.
- `stable` — actively maintained, agent can rely on it.
- `deprecated` — superseded by code change; do not use.

Do not invoke an agent against a `deprecated` skill.

### Snapshot regeneration

When a builder change legitimately changes output:

```bash
poetry run pytest hvantk/tests/test_<source>_builder.py --regenerate-snapshots
# For Hail-marked builders only:
poetry run pytest hvantk/tests/test_<source>_builder.py -m hail --regenerate-snapshots
```

Review the diff in `hvantk/tests/snapshots/<source>/` before committing. The regenerated snapshot is the new ground truth.

### Skill-driven changes

Skill-assisted code changes (where an agent generates the builder) are reviewed identically to hand-written code. NEVER auto-merge a skill-generated PR. Always run the full test suite, especially the affected `test_<source>_builder.py`.
