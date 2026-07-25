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
    clean: tuple[str, ...]
    conflicted: dict[str, frozenset[str]]
    unknown: tuple[str, ...]

    @property
    def all_columns(self) -> tuple[str, ...]:
        return tuple(self.clean) + tuple(self.conflicted) + tuple(self.unknown)


def _classes(sources, equivalence) -> frozenset[str]:
    """Map raw source names onto their equivalence classes.

    A source with no declared class maps to itself, so an unrecognised name still
    conflicts with an identical name and with nothing else.
    """
    member_to_class = {
        member: klass for klass, members in equivalence.items() for member in members
    }
    return frozenset(member_to_class.get(s, s) for s in sources)


def resolve_arms(feature_provenance, label_provenance, equivalence) -> ArmAssignment:
    """Split feature columns into clean / conflicted / unknown for a given label source.

    ``feature_provenance`` maps column -> the set of sources the predictor was trained on,
    or ``None`` when undeclared. An empty set means "trained on nothing label-derived" and
    is clean; ``None`` means "nobody said", which is NOT the same thing and must not be
    treated as clean -- a plugin author who declares nothing gets the conservative default.
    """
    label_classes = _classes(label_provenance, equivalence)
    clean, conflicted, unknown = [], {}, []
    for column, sources in feature_provenance.items():
        if sources is None:
            unknown.append(column)
            continue
        overlap = _classes(sources, equivalence) & label_classes
        if overlap:
            conflicted[column] = overlap
        else:
            clean.append(column)
    return ArmAssignment(tuple(clean), conflicted, tuple(unknown))
