# Synthetic UCSC Cell Browser `adult-ctx` fixture

SYNTHETIC. Contains no real expression values or cell counts from the "Adult
Cortex Meta-Atlas" collection (`adult-ctx-meta-atlas` on https://cells.ucsc.edu/).

## Why synthetic

A review round (#426) found:

- `https://genome.ucsc.edu/license/`, previously cited for this fixture, covers
  "the data files and database tables used by the Genome Browser" — it does
  not cover UCSC Cell Browser datasets, which are a separate product.
- The collection's own `dataset.json` (fetched from
  `https://cells.ucsc.edu/adult-ctx-meta-atlas/dataset.json`) carries no
  `license` or `rights` field at all.
- The collection pools 16 distinct third-party studies (its `Dataset` metadata
  field has 16 distinct values), each potentially under its own terms that
  UCSC does not restate.

An absent licence limits what a fixture may contain, not whether one exists
(`_conventions` SKILL.md §9), so the real pseudo-bulk summary previously
committed here was replaced with a synthetic one.

## What is kept vs. fabricated

- **Kept exactly**: the 19 `Class` group labels (column headers in
  `expression_matrix-synthetic.tsv` / row labels in `metadata-synthetic.tsv`:
  `Astrocyte`, `Endothelial`, `IT`, `L4 IT`, `L5 ET`, `L5/6 IT Car3`, `L5/6 NP`,
  `L6 CT`, `L6b`, `LAMP5`, `Microglia`, `OPC`, `Oligodendrocyte`, `PAX6`,
  `PVALB`, `Pericyte`, `SST`, `VIP`, `VLMC`) and the 250 gene identifiers (row
  labels in `expression_matrix-synthetic.tsv`) — these are referenced by the
  drift-probe docs in `hvantk/skills/ucsc_cellbrowser/SKILL.md` §12 and are
  plain identifiers, not derived data.
- **Fabricated**: every expression value in `expression_matrix-synthetic.tsv`
  and every `n_cells` count in `metadata-synthetic.tsv`.

## Generation

Seeded RNG (seed `20261013`), styled like the real file's value conventions
(exact zeros for most entries, small decimals, occasional scientific
notation). Regenerated `../../snapshots/ucsc-cellbrowser-adult-ctx/*.json` via
`pytest hvantk/skills/ucsc_cellbrowser/tests/test_builder.py -k adult-ctx --regenerate-snapshots`.
