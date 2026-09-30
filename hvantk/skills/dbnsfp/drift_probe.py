"""dbnsfp drift probe: release-list scrape of the dbnsfp.org releases page.

Follow-up to issue #321 and probe_version 1/2: this probe used to watch the legacy
Google Sites landing page (``https://sites.google.com/site/jpopgen/dbNSFP``), which
now tells readers to access the current dbNSFP at dbNSFP.org and is frozen at v4.9
(2024-08-08, verified 2026-09-30 -- the page still advertises only up to
``dbNSFPv4.9*.zip``). dbNSFP has since shipped v5.1 through v5.4 and the probe
reported clean throughout, because the page it watched could no longer change.
This version points at ``https://www.dbnsfp.org/releases/`` instead, which lists
every release as it happens.

**Why the release list stays under ``headers`` (a schema signal) rather than moving
to ``extras`` (routine).** A new dbNSFP release is, from the releases page's point of
view, routine upstream activity -- exactly the kind of signal `_conventions` section
12 says belongs in `extras`. But that guidance is conditional: *"A pure content signal
belongs in extras, but only once the probe can actually see the schema"* -- moving it
while the probe stays schema-blind converts an uninformative positive into a silent
false negative (issue #271). This probe is permanently schema-blind: the real TSV's
column header lives inside the ~50 GB academic release, obtained only through
institutional-email registration (SKILL.md section 2), so there is no way to fetch it
on a cron. Contrast clinvar's probe (issue #333), which moved its own version-bump
signal to `extras` -- but only after rewriting itself to parse the live VCF's actual
``##INFO``/``##FORMAT`` header as a true schema signal in `headers`/`checksums`, so a
real schema change is still caught independently of the routine content bump. dbNSFP
has no equivalent fallback, so the release list is the only signal available and it
stays a schema signal: every new release keeps opening its own ``drift:schema`` PR
rather than being silently batched into ``drift:routine``.

**The raw page body must still never be hashed.** A same-day check (2026-09-30) found
two consecutive fetches byte-identical, unlike the old Google Sites page. But the
footer carries a copyright year range (``2024–2026`` as of this writing) that
rolls forward every January with no release involved -- the same rendering-noise
failure mode as Google Sites, just on a yearly cadence instead of a per-request one.
The extracted release list is compared instead, exactly as before.

**Release naming on the new page.** Most releases print as a bare version shared by
both branches, e.g. ``dbNSFP v5.4 (August 1, 2026)``: one academic build, one
commercial build, same version number. Occasionally a patch ships commercial-only,
printed with the branch suffix glued directly onto the version, e.g.
``dbNSFP v5.1.1c (April 24, 2025)`` -- confirmed live on 2026-09-30, and the entry
carries no companion ``5.1.1a``. Such an entry must still appear in the release list
(a real event worth a schema PR) but must not be reported as the latest *academic*
release, since no academic file shipped that round.

What this detects: a new dbNSFP release being advertised on the releases page,
academic or commercial, including a commercial-only patch. What it cannot detect: an
in-place change to a release's actual file contents or column schema -- the page
names the release either way, and the archive itself is never fetched. Those limits
are recorded in SKILL.md.
"""

from __future__ import annotations

import re
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError
from hvantk.core.utils.http import request_with_retry

PROBE_VERSION = 3
DBNSFP_RELEASES_URL = "https://www.dbnsfp.org/releases/"

# Every release -- current or past -- renders as one of these spans (verified live
# 2026-09-30); the probe does not need to tell the "Current Release" <h1> section
# apart from the "Past Releases" <h2> one, since every entry carries this same class
# regardless of which section it sits in.
_RELEASE_TITLE_REGEX = re.compile(r'<span class="release-title">([^<]+)</span>')

# "dbNSFP v5.4 (August 1, 2026)" -> version="5.4", date="August 1, 2026".
# "dbNSFP v5.1.1c (April 24, 2025)" -> version="5.1.1c" (branch suffix glued on).
# The non-greedy version group still consumes through a trailing letter suffix,
# because nothing shorter lets the rest of the pattern match -- there is no
# whitespace between the digits and the suffix for it to stop at early.
_VERSION_DATE_REGEX = re.compile(
    r"^dbNSFP\s+v([0-9][\w.]*?)\s*\(([^)]+)\)$", re.IGNORECASE
)

# An explicit trailing "c" means commercial-only: confirmed live as the shape of the
# v5.1.1c patch, which has no accompanying academic build. Excluded from "latest
# academic" below regardless of how high its version number is.
_COMMERCIAL_ONLY_REGEX = re.compile(r"^(\d+(?:\.\d+)*)c$", re.IGNORECASE)

# A bare version (the common case: one build, both branches) or an explicit trailing
# "a" both mean an academic file shipped that round. Unlike the legacy landing page,
# where every release always listed distinct dbNSFPv<x>a.zip / dbNSFPv<x>c.zip
# archives and the suffix was mandatory, the new page's shared releases carry no
# suffix at all -- so, unlike probe_version 2, a bare version must count as academic.
_ACADEMIC_TOKEN_REGEX = re.compile(r"^(\d+(?:\.\d+)*)a?$", re.IGNORECASE)

_FILENAME = "dbNSFP-release-index"
_TIMEOUT_S = (5.0, 15.0)


def _numeric_parts(version: str) -> tuple[int, ...]:
    return tuple(int(p) for p in version.split("."))


def _latest_academic(versions: list[str]) -> str | None:
    """Highest version with an academic build, compared componentwise.

    ``versions`` are the raw version tokens exactly as printed on the page (e.g.
    "5.4", "5.1.1c"). A commercial-only token is skipped outright, never entering the
    comparison regardless of its numeric value -- a higher-numbered commercial-only
    patch must not shadow the last real academic release.
    """
    eligible = []
    for v in versions:
        if _COMMERCIAL_ONLY_REGEX.match(v):
            continue
        m = _ACADEMIC_TOKEN_REGEX.match(v)
        if m:
            eligible.append((_numeric_parts(m.group(1)), m.group(1)))
    if not eligible:
        return None
    # Componentwise so 5.10 beats 5.4; a text sort would not.
    return max(eligible)[1]


def fetch_fingerprint() -> dict:
    """Fingerprint the release entries advertised on the dbNSFP releases page."""
    try:
        resp = request_with_retry(
            "GET", DBNSFP_RELEASES_URL, timeout=_TIMEOUT_S, allow_redirects=True
        )
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    # Decoded explicitly rather than via `resp.text`, for the same reason as before:
    # requests' charset fallback for text/html is environment dependent when no
    # charset is declared. This page does declare UTF-8, but there is no upside to
    # relying on that where an explicit decode costs nothing.
    body = resp.content.decode("utf-8", errors="replace")

    titles = [t.strip() for t in _RELEASE_TITLE_REGEX.findall(body)]
    # Fail closed. No release-title spans means the page moved or its markup
    # changed, not that dbNSFP shipped zero releases; recording it would bake an
    # empty baseline that every later run compares equal to.
    if not titles:
        raise DriftProbeError(
            "dbNSFP releases page contained no release-title entries; "
            "the page layout has probably changed."
        )

    versions: list[str] = []
    for title in titles:
        m = _VERSION_DATE_REGEX.match(title)
        # Also fail closed. A title that does not fit "dbNSFP v<version> (<date>)"
        # means the page's release format changed; silently dropping it would shrink
        # the projection instead of reporting the real shape.
        if not m:
            raise DriftProbeError(
                f"dbNSFP release entry {title!r} did not match the expected "
                "'dbNSFP v<version> (<date>)' shape; the page layout has probably "
                "changed."
            )
        versions.append(m.group(1))

    latest = _latest_academic(versions)
    # Also fail closed here. Silently recording source_version: null would let the
    # bot commit that null as the baseline, after which the probe reports clean
    # forever having quietly stopped identifying a release at all.
    if latest is None:
        raise DriftProbeError(
            f"dbNSFP advertised {len(versions)} release(s) but none carried an "
            "academic build (all were explicit commercial-only 'c' patches); the "
            "release scheme has probably changed."
        )

    return {
        "probe_version": PROBE_VERSION,
        "source_version": latest,
        # The release list IS the signal -- see the module docstring for why it stays
        # under `headers` (a SCHEMA_KEY to `.github/scripts/drift_to_pr.py`) rather
        # than `extras`. A sha256 over it would be a pure function of a value already
        # in the compared surface.
        "headers": {_FILENAME: sorted(set(titles))},
        "checksums": {},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
