"""Leakage tiers for the columns of an ML training dataset.

A tier says how much a column can leak the label or allele frequency into a model:

- **L1** -- trained on nothing label-derived (conservation, structure, expression, ...).
- **L2** -- trained on the same curated disease labels the dataset is labelled with
  (ClinVar, HGMD, OMIM, ...), or declared so by the spec (e.g. predictors trained on
  simulated vs observed alleles, which are AF-biased).
- **L3** -- trained on population allele frequency, or a gene/region constraint score
  built from population depletion.
- **L4** -- a population allele frequency, or a score that takes one as an input.

Nested benchmark arms are L1, <= L2, <= L3 and all. Tiers are derived, not hand-assigned:
from each column's declared ``trained_on`` (plugin manifests) and the spec's explicit
declarations, first matching rule wins:

1. declared ``L4`` -> L4 (``trained_on`` describes training data, not inputs);
2. declared ``L2``, or ``trained_on`` conflicts with the label's sources -> L2;
3. declared ``L3``, or ``trained_on`` contains ``population_af`` -> L3;
4. declared ``L1``, or any declared ``trained_on`` -> L1;
5. otherwise the column is undeclared: excluded (``None``).

A declaration can therefore never hide circularity: a column trained on ClinVar is L2
even when the spec says L1. Pure Python, no Hail.
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
    # here must mean "undeclared" (rule 5), not L2.
    declared_sources = {
        item.column: item.trained_on for item in inputs if item.trained_on is not None
    }
    conflicted = set(
        resolve_arms(declared_sources, frozenset(label_sources), equivalence).conflicted
    )

    tiers: dict[str, str | None] = {}
    for item in inputs:
        sources = item.trained_on
        if item.declared == "L4":
            tier = "L4"
        elif item.declared == "L2" or item.column in conflicted:
            tier = "L2"
        elif item.declared == "L3" or (sources is not None and POPULATION_AF in sources):
            tier = "L3"
        elif item.declared == "L1" or sources is not None:
            tier = "L1"
        else:
            tier = None
        tiers[item.column] = tier
    return tiers
