"""Resolve a no-parse dataset's builder input to the single raw file it needs.

A plugin dataset that declares no ``lifecycle.parse`` has nothing to hand its builder
but the raw download directory: ``hvantk reprocess`` forwards ``raw_dir`` verbatim as
``parsed_input`` in that case ("No parse stage declared: the builder consumes the raw
dir directly", ``hvantk/tools/plugins/reprocess_cli.py``). A builder that reads a
single file -- ``hl.import_table``, ``hl.import_vcf``, plain ``open()`` -- cannot
consume a directory, so it must resolve ``parsed_input`` to the one file inside it
before doing anything else.
"""

from __future__ import annotations

import datetime
import fnmatch
import logging
import os
import re

logger = logging.getLogger(__name__)

# An ISO 8601 calendar date (YYYY-MM-DD) in a filename, e.g. the `<prefix>-2026-10-14.csv`
# clingen and gencc each download to. The lookaround on both sides requires the match not
# be embedded in a longer digit run (so a stray 9- or 11-digit number never looks like a
# date), and `_dated_name_or_none` below additionally validates the matched substring is a
# real calendar date, not just a string shaped like one (e.g. `2026-13-45`).
_ISO_DATE_RE = re.compile(r"(?<!\d)\d{4}-\d{2}-\d{2}(?!\d)")


def _dated_name_or_none(name: str) -> tuple[str, str] | None:
    """Split *name* into ``(stem_with_date_removed, date_string)``, or None.

    None when *name* carries no ``_ISO_DATE_RE`` match, or when it does but the
    matched substring is not a real calendar date (``datetime.date.fromisoformat``
    raises) -- either way, *name* does not count as "dated" for `latest_dated`
    below. ``stem_with_date_removed`` lets the caller check that every candidate is
    the *same* file published on a different day, not merely that it happens to
    carry *some* date: an unrelated dated file (a different report that happens to
    also be named `<other>-<date>.tsv`) must not take part.
    """
    match = _ISO_DATE_RE.search(name)
    if match is None:
        return None
    date_str = match.group(0)
    try:
        datetime.date.fromisoformat(date_str)
    except ValueError:
        return None
    stem = name[: match.start()] + name[match.end() :]
    return stem, date_str


def resolve_raw_file(
    parsed_input, *, patterns: tuple[str, ...] = ("*",), latest_dated: bool = False
) -> str:
    """Resolve *parsed_input* to a single file path.

    A path that is not a directory is returned unchanged (as ``str``), which keeps
    direct ``build(<file>)`` calls (tests, snapshots) working exactly as before.

    A directory is resolved to the single non-hidden regular file inside it whose
    basename matches one of *patterns* (``fnmatch``, case-insensitive, applied to the
    immediate children only -- it does not recurse). Pass patterns that match both
    what the dataset's downloader actually writes and what its test fixture is named,
    so the same call site resolves a real download directory and the committed
    fixture directory alike.

    Parameters
    ----------
    parsed_input : str | os.PathLike
        The builder's ``parsed_input``: a file path, or -- for a no-parse dataset
        reached through ``hvantk reprocess`` -- the raw download directory.
    patterns : tuple[str, ...]
        ``fnmatch`` patterns the target file's basename must match (matched
        case-insensitively). Defaults to ``("*",)`` (any non-hidden regular file),
        which is correct for a dataset whose raw dir only ever holds one file
        regardless of its name.
    latest_dated : bool
        When several files match *patterns* and this is True, resolve to the one
        whose name carries the latest real ISO 8601 date (``YYYY-MM-DD``,
        validated with ``datetime.date.fromisoformat`` -- a date-shaped but
        invalid substring such as ``2026-13-45`` does not count) -- but only when
        ALL of these also hold:

        * every match carries one (a mix of dated and undated names, e.g.
          clinvar's undated ``clinvar.vcf.gz`` next to a dated archive, still
          raises: there is then no uniform rule for picking one);
        * every match is the name of the *same* file with its date removed (an
          unrelated dated file that happens to match *patterns* too -- a
          different report, also named ``<other>-<date>.tsv`` -- must not take
          part, or it could silently be picked over, or in place of, the real
          series);
        * the latest date is not itself a tie across more than one candidate.

        Any of those failing falls through to the same "found N" error as when
        *latest_dated* is False. Set this for a dataset whose downloader labels
        each run's output with today's date (clingen, gencc), so a second day's
        `reprocess` run -- which leaves yesterday's file in place alongside
        today's -- does not turn "found 2" into a permanent build failure. The
        ignored older file(s) are logged at WARNING, not silently dropped.

    Returns
    -------
    str
        The resolved file path.

    Raises
    ------
    ValueError
        If the directory holds zero, or more than one (and *latest_dated* does not
        apply), file matching *patterns*. The message names the directory, the
        patterns, and the files found, and says what to do.
    """
    path = str(parsed_input)
    if not os.path.isdir(path):
        return path

    candidates = sorted(
        name
        for name in os.listdir(path)
        if not name.startswith(".")
        and os.path.isfile(os.path.join(path, name))
        and any(
            fnmatch.fnmatchcase(name.lower(), pattern.lower()) for pattern in patterns
        )
    )

    if len(candidates) == 1:
        return os.path.join(path, candidates[0])

    if latest_dated and len(candidates) > 1:
        dated = [(name, _dated_name_or_none(name)) for name in candidates]
        if all(info is not None for _, info in dated) and (
            len({info[0] for _, info in dated}) == 1
        ):
            best_date = max(info[1] for _, info in dated)
            winners = [name for name, info in dated if info[1] == best_date]
            if len(winners) == 1:
                chosen_name = winners[0]
                ignored = sorted(name for name in candidates if name != chosen_name)
                logger.warning(
                    "resolve_raw_file: %d files matching %s in %r are the same "
                    "file under different dates; using the latest (%r, dated %s) "
                    "and ignoring %s",
                    len(candidates),
                    patterns,
                    path,
                    chosen_name,
                    best_date,
                    ignored,
                )
                return os.path.join(path, chosen_name)

    raise ValueError(
        f"expected exactly one file matching {patterns} in directory {path!r}, "
        f"found {len(candidates)}: {candidates}. Leave exactly one matching file in "
        "--raw-dir, or build from the file directly: `hvantk reprocess ... "
        "--skip-download --skip-parse --intermediate <file>` (--raw-dir itself must "
        "be a directory), or, from the Python API, pass the file path as parsed_input."
    )
