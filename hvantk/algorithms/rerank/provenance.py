"""Decide which feature columns may contribute to the headline result.

Circularity is a property of a feature/label PAIR, not of a feature. REVEL (trained on
HGMD/ClinVar) is badly circular against a ClinGen- or GenCC-derived label, and far less so
against a purely burden-derived one. A static per-feature "leakage tier" cannot express
that, so provenance is declared on both sides and the conflict is computed.

This matters because statistical feature selection cannot detect circularity -- it REWARDS
it. A univariate filter ranks REVEL highly precisely because REVEL was trained on genes
like these. Provenance is therefore not recoverable from the data and has to be declared.
"""
from __future__ import annotations

from dataclasses import dataclass

# Sources within one class conflict with each other; sources in different classes never do.
# Extending this is a data edit in selection.yaml, not a code change.
DEFAULT_EQUIVALENCE: dict[str, list[str]] = {
    "curated_disease_db": ["ClinVar", "ClinGen", "GenCC", "HGMD", "OMIM"],
}


@dataclass(frozen=True)
class ArmAssignment:
    """Which columns may carry the headline result, for one feature/label pair.

    ``conflicted`` maps a barred column to the equivalence classes it conflicts on.
    ``undeclared`` is a *subset of* ``conflicted``'s keys: columns barred not because a
    conflict was computed but because nobody declared what the predictor was trained on.
    They are kept nameable so the manifest gap can be closed; without the list, an
    undeclared column is indistinguishable from a genuinely circular one and no one can
    tell which declaration is missing.
    """

    clean: tuple[str, ...]
    conflicted: dict[str, frozenset[str]]
    undeclared: tuple[str, ...]

    @property
    def all_columns(self) -> tuple[str, ...]:
        # `undeclared` is not added: its members are already keys of `conflicted`.
        return tuple(self.clean) + tuple(self.conflicted)


def _classes(sources, equivalence) -> frozenset[str]:
    """Map raw source names onto their equivalence classes.

    A source with no declared class maps to itself, so an unrecognised name still
    conflicts with an identical name and with nothing else.

    A source listed under TWO classes is rejected rather than resolved. Building the
    lookup as a plain dict comprehension would let the last class processed win by
    iteration order, silently reclassifying that source and moving columns between the
    clean and all arms -- a wrong headline, not an error. Same rule the CLI already
    applies to a repeated ``--prepared`` axis (PR #229): reject the ambiguity, do not
    guess at it.
    """
    member_to_class: dict[str, str] = {}
    for klass, members in equivalence.items():
        for member in members:
            prior = member_to_class.get(member)
            if prior is not None and prior != klass:
                raise ValueError(
                    f"equivalence source {member!r} is declared in two classes "
                    f"({prior!r} and {klass!r}); a source belongs to exactly one class, "
                    "because its class decides which arm a feature lands in. Remove the "
                    "duplicate."
                )
            member_to_class[member] = klass
    return frozenset(member_to_class.get(s, s) for s in sources)


def resolve_arms(feature_provenance, label_provenance, equivalence) -> ArmAssignment:
    """Split feature columns into clean / conflicted for a given label source.

    ``feature_provenance`` maps column -> the set of sources the predictor was trained on,
    or ``None`` when undeclared. An empty set means "trained on nothing label-derived" and
    is clean; ``None`` means "nobody said", which is NOT the same thing and must not be
    treated as clean -- a plugin author who declares nothing gets the conservative default.

    That conservative default is **conflicted**, presumed to conflict on whatever classes
    the label itself derives from: if we do not know what a predictor was trained on, the
    assumption that costs us signal is safer than the one that inflates the headline.
    Undeclared columns were previously reported in a third bucket of their own, which
    barred them from the clean arm just the same but described them as merely
    unclassified -- a column silently treated as circular ought to say so. Their names are
    still carried, in ``undeclared``, so the missing declaration can be found and written.

    Arm membership is unchanged by this: the caller builds the clean arm from
    ``assignment.clean`` and leaves ``all`` unrestricted, so a column that moved from the
    old ``unknown`` bucket into ``conflicted`` was excluded from clean and present in all
    both before and after. No result computed under the previous behaviour shifts.
    """
    label_classes = _classes(label_provenance, equivalence)
    clean, conflicted, undeclared = [], {}, []
    for column, sources in feature_provenance.items():
        if sources is None:
            conflicted[column] = label_classes
            undeclared.append(column)
            continue
        overlap = _classes(sources, equivalence) & label_classes
        if overlap:
            conflicted[column] = overlap
        else:
            clean.append(column)
    return ArmAssignment(tuple(clean), conflicted, tuple(undeclared))
