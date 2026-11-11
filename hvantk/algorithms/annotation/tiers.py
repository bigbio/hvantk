"""Leakage tiers for the columns of an ML training dataset.

A tier says how much a column can leak the label or allele frequency into a model:

- **L1** -- trained on nothing label-derived (conservation, structure, expression, ...).
- **L2** -- trained on the same curated disease labels the dataset is labelled with
  (ClinVar, HGMD, OMIM, ...), or declared so by the spec (e.g. predictors trained on
  simulated vs observed alleles, which are AF-biased).
- **L3** -- trained on population allele frequency, or a gene/region constraint score
  built from population depletion.
- **L4** -- a population allele frequency, or a score that takes one as an input.

Nested benchmark arms are L1, <= L2, <= L3 and all, so a higher tier keeps a column out of
more arms. Tiers are derived, not hand-assigned. From a column's declared ``trained_on``
(plugin manifests):

- ``trained_on`` conflicts with the label's sources -> at least L2;
- ``trained_on`` contains ``population_af`` -> at least L3;
- any other declared ``trained_on`` (including an empty one) -> L1.

The column's tier is the higher of that and the spec's ``tier:`` declaration, so a
declaration can raise a tier but never lower it: a column trained on ClinVar is at least L2
whatever the spec says, and a declared L4 (an allele frequency, or a score taking one as
input; ``trained_on`` describes training data, not inputs) stays L4. A column with neither
a declared ``trained_on`` nor a ``tier:`` is undeclared: excluded (``None``). Pure Python,
no Hail.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

from hvantk.algorithms.rerank.provenance import DEFAULT_EQUIVALENCE, resolve_arms

TIERS = ("L1", "L2", "L3", "L4")
POPULATION_AF = "population_af"


@dataclass(frozen=True)
class TierInput:
    """What is known about one column.

    ``trained_on`` is ``None`` when nobody declared the column's training data, and an
    empty set when it was declared to be trained on nothing label-derived.
    """

    column: str
    trained_on: frozenset[str] | None
    declared: str | None = None


def resolve_tiers(
    inputs: Sequence[TierInput],
    label_sources: Iterable[str],
    equivalence: Mapping[str, list[str]] = DEFAULT_EQUIVALENCE,
) -> dict[str, str | None]:
    """Return column -> tier (``"L1"``..``"L4"``), or ``None`` for an excluded column.

    Raises
    ------
    ValueError
        If a declared tier is not one of :data:`TIERS`.
    """
    for item in inputs:
        if item.declared is not None and item.declared not in TIERS:
            raise ValueError(
                f"column {item.column!r} declares tier {item.declared!r}; "
                f"expected one of {', '.join(TIERS)}"
            )
    # Only declared columns go to resolve_arms: it reads None as "conflicted", which
    # here must mean "undeclared" (excluded), not L2.
    declared_sources = {
        item.column: item.trained_on for item in inputs if item.trained_on is not None
    }
    conflicted = set(
        resolve_arms(declared_sources, frozenset(label_sources), equivalence).conflicted
    )

    tiers: dict[str, str | None] = {}
    for item in inputs:
        sources = item.trained_on
        candidates = [item.declared] if item.declared is not None else []
        if sources is not None:
            candidates.append("L1")
            if item.column in conflicted:
                candidates.append("L2")
            if POPULATION_AF in sources:
                candidates.append("L3")
        tiers[item.column] = max(candidates, key=TIERS.index) if candidates else None
    return tiers
