"""`normalize_etag` must reduce every real ETag form to one bare tag.

The shipped-once bug was `str.strip('"')`, which strips a character SET from both
ends: `W/"abc"` kept its `W/` prefix because the leading `W` blocked the left
strip. Because the drift bot regenerates drifted baselines automatically, one
mangled value bakes in permanently.
"""

from __future__ import annotations

import pytest

from hvantk.core.plugin.api import normalize_etag


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ('"abc123"', "abc123"),          # strong
        ('W/"abc123"', "abc123"),        # weak -- the strip() bug
        ('w/"abc123"', "abc123"),        # weak, lowercase
        ('"abc123"-gzip', "abc123"),     # transform suffix appended by mod_deflate
        ('  "abc123"  ', "abc123"),      # surrounding whitespace
        ("abc123", "abc123"),            # unquoted, seen from some proxies
    ],
)
def test_every_form_reduces_to_the_same_tag(raw, expected):
    assert normalize_etag(raw) == expected


@pytest.mark.parametrize("raw", [None, "", '""', "   ", 'W/""'])
def test_absent_or_empty_becomes_none(raw):
    """Callers fail closed on None. Returning '' instead would be recorded as an
    empty digest, which placeholder_baseline_reason later reads as a hand-seeded
    baseline -- an unbreakable probe_failed loop."""
    assert normalize_etag(raw) is None


def test_the_old_strip_bug_would_fail_these():
    """Documents precisely what regressed, so a revert to strip('"') is caught."""
    assert 'W/"abc"'.strip('"') == 'W/"abc'      # the bug
    assert normalize_etag('W/"abc"') == "abc"     # the fix
