"""Paralogue blocks for cross-validation, and the one case where blocking must refuse.

Every number the rerank engine produced used random 5-fold CV. Gene families break that:
paralogues share sequence, share constraint, share expression pattern, and share disease
status, so a family split across train and test lets the model recognise a relative rather
than generalise. Blocking is the standard answer, and nothing in the library offered it.

**Primary group, not connected components.** HGNC's ``gene_group`` is multi-membership and
pipe-separated, and the obvious blocker -- union-find over shared membership -- is not merely
conservative here, it is INVALID. Membership chains transitively (A in {X,Y}, B in {Y,Z}, C in
{Z,W}, ...), so the closure can collapse a large fraction of a gene universe into one
component. ``StratifiedGroupKFold`` must place a whole block in one fold, so the pooled
out-of-fold AUC then estimates something different from the random-fold AUC it is compared
against; the tell-tale is a "harder" grouped model scoring HIGHER than random folds, which
blocking cannot do. Blocking on the first-listed (primary) group keeps blocks small; the
residual leak -- a pair sharing only a secondary group -- is a far smaller error than an
incomparable estimator.

**The ceiling is a hard abort.** Same argument, one step on: if the primary group is itself
oversized for this universe, the run must fail rather than produce a plausible number. A
warning is not enough precisely because the failure is silent -- the output looks fine and the
tell-tale is a counterintuitive AUC that a reader has no reason to question.
"""
from __future__ import annotations

import hashlib
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

DEFAULT_MAX_BLOCK_FRAC = 0.10
"""Largest permitted block, as a fraction of the universe -- a guard rail rather than a
tuned threshold."""


class DominantBlockError(ValueError):
    """One block holds too much of the universe for blocked folds to mean anything."""


def _validate_max_block_frac(value) -> float:
    """Shared by ``BlockPolicy`` and ``gene_blocks``, which each validate independently --
    the former at policy-construction time, the latter wherever it is called directly."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"max_block_frac must be a real number, not bool; got {value!r}")
    if not math.isfinite(value):
        raise ValueError(f"max_block_frac must be a finite number; got {value!r}")
    if not 0.0 < float(value) <= 1.0:
        raise ValueError(f"max_block_frac must be in (0, 1]; got {value}")
    return float(value)


@dataclass(frozen=True)
class BlockPolicy:
    """Where the gene groups come from, and how large a block may be.

    ``table`` is REQUIRED, not defaulted: an empty policy would resolve to an empty
    gene-group mapping, every gene would become its own singleton block, and "blocked" CV
    would silently equal random folds while a null built from it still recorded a blocking
    -- a control that looks on but is off, which is exactly what this batch exists to
    prevent. It is a path rather than a packaged resource because the file is versioned
    upstream and updated independently of this library.
    """

    table: str | os.PathLike
    max_block_frac: float = DEFAULT_MAX_BLOCK_FRAC

    def __post_init__(self) -> None:
        raw = self.table
        table = os.fspath(raw) if isinstance(raw, os.PathLike) else raw
        if not isinstance(table, str) or not table.strip():
            raise ValueError(
                f"BlockPolicy.table must be a non-empty path (str or os.PathLike); got "
                f"{raw!r}"
            )
        object.__setattr__(self, "table", table)
        object.__setattr__(
            self, "max_block_frac", _validate_max_block_frac(self.max_block_frac)
        )


def block_digest(labels) -> str:
    """A content hash of block labels, in the SAME gene order they were computed in.

    Two nulls generated under different block tables -- another HGNC release, another
    grouping -- must not compare equal just because both happen to be "blocked": this is
    what lets ``ControlSetting`` tell them apart rather than silently answering for a
    blocking it was not built under.
    """
    return hashlib.sha1(np.asarray(labels, dtype=np.int64).tobytes()).hexdigest()


@dataclass(frozen=True)
class BlockReport:
    """The blocking, plus what a reader needs to judge whether it binds.

    Grouped folds cannot be balanced the way stratified random folds are, so how much of the
    matrix is actually constrained is part of the result, not a debug detail.
    """

    blocks: np.ndarray
    n_blocks: int
    largest: int
    largest_frac: float
    n_in_multi: int

    @property
    def digest(self) -> str:
        """See :func:`block_digest`."""
        return block_digest(self.blocks)


def load_gene_groups(path) -> dict:
    """``{symbol: gene_group}`` from an HGNC complete-set TSV, Approved entries only.

    Withdrawn and merged entries are dropped rather than kept: a withdrawn symbol's family
    is not a statement about the gene that carries that symbol today.
    """
    import pandas as pd

    path = Path(path)
    frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    missing = [c for c in ("symbol", "gene_group") if c not in frame.columns]
    if missing:
        raise ValueError(
            f"{path}: gene-group table is missing column(s) {missing}; expected an HGNC "
            "complete-set TSV with 'symbol' and 'gene_group' (and optionally 'status')"
        )
    if "status" in frame.columns:
        frame = frame[frame["status"] == "Approved"]
    return dict(zip(frame["symbol"], frame["gene_group"].fillna("")))


def gene_blocks(
    genes: Sequence[str],
    gene_group: Mapping,
    *,
    max_block_frac: float = DEFAULT_MAX_BLOCK_FRAC,
) -> BlockReport:
    """Block codes aligned to ``genes``, one block per primary HGNC gene group.

    Raises ``DominantBlockError`` when the largest block exceeds ``max_block_frac`` of the
    universe -- see the module docstring for why that is an abort rather than a warning.
    Raises ``ValueError`` when ``genes`` is non-empty and none of them appears in
    ``gene_group``: the table's keys probably do not match this universe's gene keys (e.g.
    Ensembl IDs vs HGNC symbols), and blocked folds would otherwise silently be unblocked
    with nothing to show for it. A gene present with an empty group still counts as matched.
    """
    max_block_frac = _validate_max_block_frac(max_block_frac)

    genes = list(genes)
    if genes and not any(g in gene_group for g in genes):
        raise ValueError(
            f"none of the {len(genes)} gene(s) in this universe appear in the gene-group "
            f"table (e.g. {genes[:3]!r}); the table's keys probably do not match this "
            "universe's gene keys (e.g. Ensembl IDs vs HGNC symbols), and blocked folds "
            "would silently be unblocked"
        )

    codes: dict = {}
    out = []
    for i, gene in enumerate(genes):
        raw = gene_group.get(gene)
        primary = raw.split("|")[0].strip() if isinstance(raw, str) else ""
        # An ungrouped gene has no paralogue to leak through, so it is its own block --
        # unconstrained, which is the correct treatment, not pooled with the other orphans.
        key = f"fam:{primary}" if primary else f"solo:{i}"
        out.append(codes.setdefault(key, len(codes)))

    blocks = np.asarray(out, dtype=int)
    sizes = np.bincount(blocks) if len(blocks) else np.zeros(0, dtype=int)
    largest = int(sizes.max()) if sizes.size else 0
    n = len(genes)
    if n and largest > max_block_frac * n:
        raise DominantBlockError(
            f"the largest paralogue block holds {largest}/{n} units "
            f"({largest / n:.1%}), over the {max_block_frac:.0%} ceiling "
            f"(max_block_frac). StratifiedGroupKFold must place a whole block in one fold, "
            f"so one fold would hold that entire block and the pooled out-of-fold AUC would "
            f"not estimate the same quantity as the unblocked run. Raise max_block_frac "
            f"deliberately if you accept that, or restrict the universe."
        )
    n_in_multi = int((sizes[blocks] > 1).sum()) if sizes.size else 0
    return BlockReport(
        blocks=blocks,
        n_blocks=int(sizes.size),
        largest=largest,
        largest_frac=float(largest / n) if n else 0.0,
        n_in_multi=n_in_multi,
    )
