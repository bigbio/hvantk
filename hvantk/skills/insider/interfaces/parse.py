"""Pure-Python parse of the INSIDER ``H_sapiens_interfacesALL.txt`` interface table.

No Hail, no pandas-Hail bridge -- so the whole reduction is fast-testable. The builder
wraps the returned frame into a Hail Table.

Raw format (tab-separated, one row per interacting protein PAIR)::

    P1          P2          Source  P1_IRES  P2_IRES
    A0A024R087  Q92560      ECLAIR  [265]    [51,150,655,659,661,663,666]

``*_IRES`` are the interface residue positions on that side of the pair, as a bracketed
comma list that may contain ``lo-hi`` ranges (``[29,33,36-38,41]``) or be empty (``[]``).
``Source`` is ``ECLAIR`` (predicted) or ``PDB``/``I3D`` (experimental structure).

The reduction is per PROTEIN, not per pair: each row contributes to both of its
proteins, so the table is scanned once and both sides accumulated.
"""
from __future__ import annotations

import logging
from collections import defaultdict

logger = logging.getLogger(__name__)

# Sources backed by an observed structure, as opposed to ECLAIR's predictions. Kept as a
# named set because the predicted/experimental split is the axis's study-bias control.
EXPERIMENTAL_SOURCES = frozenset({"PDB", "I3D"})

_REQUIRED_COLUMNS = ("P1", "P2", "Source", "P1_IRES", "P2_IRES")


def expand_ires(field: str) -> set[int]:
    """Expand a bracketed IRES list into a set of residue positions.

    ``"[29,33,36-38]"`` -> ``{29, 33, 36, 37, 38}``; ``"[]"`` and ``""`` -> ``set()``.
    Unparseable tokens are skipped rather than raising: the file is a third-party dump
    and one malformed residue must not abort a 123k-row build.
    """
    if not field:
        return set()
    body = field.strip().strip("[]").strip()
    if not body:
        return set()
    out: set[int] = set()
    for token in body.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token[1:]:  # [1:] so a leading '-' is not read as a range separator
            lo, _, hi = token.partition("-")
            try:
                out.update(range(int(lo), int(hi) + 1))
            except ValueError:
                continue
        else:
            try:
                out.add(int(token))
            except ValueError:
                continue
    return out


def parse_interfaces(path: str):
    """Collapse the pair-level interface table to one row per protein.

    Returns a pandas DataFrame with columns ``uniprot_id``, ``n_partners``,
    ``n_partners_experimental``, ``n_partners_predicted``, ``n_interface_residues``,
    sorted by ``uniprot_id``.

    ``n_interface_residues`` is the size of the UNION of that protein's interface
    residues across all of its interactions -- a residue at two different interfaces is
    counted once, so the value is bounded by protein length and reads as "how much of
    this protein is interface".
    """
    import pandas as pd

    partners: dict[str, set[str]] = defaultdict(set)
    partners_exp: dict[str, set[str]] = defaultdict(set)
    residues: dict[str, set[int]] = defaultdict(set)

    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        missing = [c for c in _REQUIRED_COLUMNS if c not in header]
        if missing:
            raise ValueError(
                f"{path}: interface table missing column(s) {missing}; got {header}"
            )
        idx = {c: header.index(c) for c in _REQUIRED_COLUMNS}

        n_rows = 0
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            p1, p2 = f[idx["P1"]].strip(), f[idx["P2"]].strip()
            if not p1 or not p2:
                continue
            experimental = f[idx["Source"]].strip() in EXPERIMENTAL_SOURCES
            for me, other, ires_col in (
                (p1, p2, "P1_IRES"),
                (p2, p1, "P2_IRES"),
            ):
                partners[me].add(other)
                if experimental:
                    partners_exp[me].add(other)
                residues[me] |= expand_ires(f[idx[ires_col]])
            n_rows += 1

    rows = [
        {
            "uniprot_id": acc,
            "n_partners": len(partners[acc]),
            "n_partners_experimental": len(partners_exp[acc]),
            "n_partners_predicted": len(partners[acc] - partners_exp[acc]),
            "n_interface_residues": len(residues[acc]),
        }
        for acc in sorted(partners)
    ]
    logger.info("parse_interfaces: %d pair rows -> %d proteins", n_rows, len(rows))
    return pd.DataFrame(rows)
