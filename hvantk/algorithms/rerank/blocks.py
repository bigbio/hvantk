"""Paralogue blocks for cross-validation, and the one case where blocking must refuse.

Every number the rerank engine produced used random 5-fold CV. Gene families break that:
paralogues share sequence, share constraint, share expression pattern, and share disease
status, so a family split across train and test lets the model recognise a relative rather
than generalise. Blocking is the standard answer, and nothing in the library offered it.

**First-listed group, not connected components.** HGNC's ``gene_group`` is multi-membership
and pipe-separated, and the obvious blocker -- union-find over shared membership -- is not
merely conservative here, it is INVALID. Membership chains transitively (A in {X,Y}, B in
{Y,Z}, C in {Z,W}, ...), so the closure can collapse a large fraction of a gene universe into
one component. ``StratifiedGroupKFold`` must place a whole block in one fold, so the pooled
out-of-fold AUC then estimates something different from the random-fold AUC it is compared
against; the tell-tale is a blocked AUC MATERIALLY higher than the random-fold one, beyond
the across-seed spread ``--seed-sweep`` reports -- blocking is harder only ON AVERAGE, so a
single correctly blocked run can still score above random folds by chance. Blocking on the
first-listed group (in HGNC's list order -- HGNC does not designate a "primary" one) keeps
blocks small; the residual leak -- a pair sharing only a later-listed group -- is a far
smaller error than an incomparable estimator.

**The ceiling is a hard abort.** Same argument, one step on: if the first-listed group is
itself oversized for this universe, the run must fail rather than produce a plausible number.
A warning is not enough precisely because the failure is silent -- the output looks fine and
the tell-tale is a counterintuitive AUC that a reader has no reason to question.
"""

from __future__ import annotations

import hashlib
import logging
import math
import numbers
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

logger = logging.getLogger(__name__)

DEFAULT_MAX_BLOCK_FRAC = 0.10
"""Largest permitted block, as a fraction of the universe -- a guard rail rather than a
tuned threshold."""


class DominantBlockError(ValueError):
    """One block holds too much of the universe for blocked folds to mean anything."""


def _validate_max_block_frac(value) -> float:
    """Shared by ``BlockPolicy`` and ``gene_blocks``, which each validate independently --
    the former at policy-construction time, the latter wherever it is called directly.

    Accepts any ``numbers.Real`` (a bare ``float``/``int``, or a numpy scalar such as
    ``np.float32`` -- none of these are ``isinstance``-compatible with plain ``float``),
    but never a ``bool``: ``bool`` is itself a ``numbers.Real`` subtype in Python, and
    ``True``/``False`` silently reading as ``1.0``/``0.0`` is not a ceiling anyone meant to
    set. A non-number (e.g. a string) gets its own message that does not mention bool --
    conflating "wrong type entirely" with "the one wrong-type value we specifically guard
    against" would mislead a caller who passed neither.
    """
    if isinstance(value, bool):
        raise ValueError(
            f"max_block_frac must be a real number, not bool; got {value!r}"
        )
    if not isinstance(value, numbers.Real):
        raise ValueError(f"max_block_frac must be a real number; got {value!r}")
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
    -- a control that looks on but is off, which is exactly the failure mode blocking
    exists to prevent. It is a path rather than a packaged resource because the file is
    versioned upstream and updated independently of this library.
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

    Canonicalised with ``pandas.factorize`` before hashing, so the digest is a function of
    the PARTITION (which positions share a block) rather than of the raw label values:
    float or string labels work, and two labels that are merely close in value (e.g.
    ``1.2`` and ``1.7``) are never accidentally collapsed into one block the way a naive
    integer cast would collapse them (both truncate to ``1``). Hashed as fixed-width
    little-endian int64 bytes (``"<i8"``) rather than the platform's native byte order, so
    the digest does not depend on the machine it was computed on. ``gene_blocks`` already
    emits canonical 0..B-1 codes in first-appearance order, so ``factorize`` reproduces
    them unchanged and every digest computed from an engine-built blocking is unchanged
    by this.

    ``use_na_sentinel=False``, matching :func:`~hvantk.algorithms.rerank.nulls.
    _block_structure`'s convention: a NaN label is one ordinary block rather than an
    excluded sentinel (-1) in both places, so the two never disagree about how many blocks
    a NaN-containing label array has. Every label ``gene_blocks`` produces is a NaN-free
    int code, so this leaves the digest of any engine-built blocking unchanged.
    """
    import pandas as pd

    codes = pd.factorize(np.asarray(labels), use_na_sentinel=False)[0]
    return hashlib.sha1(np.asarray(codes, dtype="<i8").tobytes()).hexdigest()


@dataclass(frozen=True, eq=False)
class BlockReport:
    """The blocking, plus what a reader needs to judge whether it binds.

    Grouped folds cannot be balanced the way stratified random folds are, so how much of the
    matrix is actually constrained is part of the result, not a debug detail.

    ``eq=False``: ``blocks`` is a numpy array, whose ``==`` returns an array rather than a
    bool, which makes the dataclass-generated ``__eq__``/``__hash__`` raise the moment
    either is used -- there is no valid default to fall back to, so it is turned off rather
    than left to fail at a random call site.
    """

    blocks: np.ndarray
    n_blocks: int
    largest: int
    largest_frac: float
    n_in_multi: int
    n_matched: int
    """How many of the universe's genes were present as keys of the gene-group mapping
    (an empty group still counts as matched -- only ABSENCE from the table does not). A low
    match rate usually means the table's gene keys do not agree with the universe's (e.g. an
    older/newer HGNC release, or symbols vs. a different identifier), silently degrading
    "blocked" CV toward unblocked without any block-count signal to catch it, since a wave of
    unmatched genes still shows up as ordinary (if numerous) singletons."""

    @property
    def digest(self) -> str:
        """See :func:`block_digest`."""
        return block_digest(self.blocks)


def load_gene_groups(path) -> dict:
    """``{symbol: gene_group}`` from an HGNC complete-set TSV, Approved entries only.

    Withdrawn and merged entries are dropped rather than kept: a withdrawn symbol's family
    is not a statement about the gene that carries that symbol today. Matching against a
    universe of gene identifiers is on the approved HGNC ``symbol`` column ONLY -- aliases
    and previous symbols are not consulted, so a gene known to the universe only by an alias
    or an old symbol will not be found here and falls through to being its own singleton
    block (see ``gene_blocks``), not an error by itself.
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
    """Block codes aligned to ``genes``, one block per first-listed HGNC gene group.

    Matching a gene against ``gene_group`` is a plain key lookup -- the mapping is expected
    to already be keyed on whatever identifier ``genes`` uses (see ``load_gene_groups``,
    which keys on the approved HGNC ``symbol`` only; aliases and previous symbols are not
    matched and fall through to singletons).

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
    n_matched = 0
    for gene in genes:
        matched = gene in gene_group
        n_matched += matched
        raw = gene_group.get(gene)
        first_listed = raw.split("|")[0].strip() if isinstance(raw, str) else ""
        # An ungrouped gene has no paralogue to leak through, so it is its own block --
        # unconstrained, which is the correct treatment, not pooled with the other orphans.
        # Keyed by the gene itself (not by position): duplicate rows of the SAME gene must
        # land in the SAME singleton block, or they could be split across folds despite
        # being, trivially, each other's closest possible paralogue.
        key = f"fam:{first_listed}" if first_listed else f"solo:{gene}"
        out.append(codes.setdefault(key, len(codes)))

    n = len(genes)
    if n and n_matched < n / 2:
        logger.warning(
            "gene-group table matched only %d/%d genes (%.0f%%) in this universe; a custom "
            "family table may legitimately cover a subset, but this usually means the "
            "table's gene keys disagree with the universe's (e.g. a different HGNC release, "
            "or symbols vs. another identifier) -- the unmatched genes were treated as "
            "unconstrained singletons, not as an error",
            n_matched,
            n,
            100.0 * n_matched / n,
        )

    blocks = np.asarray(out, dtype=int)
    sizes = np.bincount(blocks) if len(blocks) else np.zeros(0, dtype=int)
    largest = int(sizes.max()) if sizes.size else 0
    # Comparing by division (largest/n > frac) rather than by multiplication
    # (largest > frac*n) makes an exact-equality boundary exact: 0.29 * 100 ==
    # 28.999999999999996 in float arithmetic, so the multiply form wrongly aborted a
    # 29/100 block at a 0.29 ceiling although the abort is pinned to strictly greater.
    # `largest > 1` first: a blocking made only of singletons is stratified random CV and
    # must never abort, but for n < 1/max_block_frac the division form alone still fires
    # (a largest of 1 exceeds a ceiling fraction of e.g. 0.02 on a 20-gene universe).
    if largest > 1 and largest / n > max_block_frac:
        code_to_key = {code: key for key, code in codes.items()}
        largest_code = int(np.argmax(sizes))
        largest_key = code_to_key[largest_code]
        # Strip whichever internal prefix built this key (see `key = ...` above) so the
        # abort message reads as a gene-group name, not an implementation-internal tag.
        # `solo:` blocks are always size 1, so `largest > 1` means this branch can never
        # actually observe one -- stripped for symmetry with `fam:`, not reachability.
        group_name = largest_key
        for prefix in ("fam:", "solo:"):
            if group_name.startswith(prefix):
                group_name = group_name[len(prefix) :]
                break
        raise DominantBlockError(
            f"the largest paralogue block (group {group_name!r}) holds {largest}/{n} units "
            f"({largest / n:.1%}), over the {max_block_frac:g} ceiling "
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
        n_matched=n_matched,
    )
