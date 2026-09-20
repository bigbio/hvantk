"""ClinVar drift probe: the VCF's own header, fetched with a ranged GET.

``probe_version`` 1 issued a bare HEAD and recorded ``Content-Length`` under ``headers``.
Because ``headers`` is a *schema* key to the drift bot (``SCHEMA_KEYS`` in
``.github/scripts/drift_to_pr.py``), every ClinVar release -- i.e. ClinVar adding
variants, its normal weekly activity -- was tiered ``drift:schema`` and given its own
pull request. The tier was therefore constant, and a constant signal carries no
information: across six regenerations in five weeks, a file that had **doubled** in size
and one that grew by 6,645 bytes produced byte-identical classifications (issue #333).

The fix is to remove the blindness rather than relabel it. Merely moving
``content_length`` into ``extras`` would have made the tiering *look* right while the
probe still could not see the schema, converting uninformative positives into genuine
false negatives -- #271's failure mode.

**How the header is obtained without the ~194 MB body.** ``clinvar.vcf.gz`` is BGZF, a
sequence of independently-decompressable gzip members, so its leading blocks decompress
on their own. A ``Range`` request for the first 64 KiB yields all 39 ``##INFO``
declarations and the ``#CHROM`` line with room to spare: measured 2026-09-21 against the
live file, the header is 45 lines / 6,406 bytes, reached within 10 whole members, and
the fetched prefix runs ~1,784 lines deep once data rows are counted. NCBI's
FTP-over-HTTPS endpoint honours ``Range`` and answers ``206``; the issue flagged this as
unverified, and it is now verified. The probe still guards against a server that ignores
it -- see ``_fetch_header_bytes``.

What moves where, and why:

* ``headers`` -- the real ``INFO``/``FORMAT`` ID lists parsed from the VCF header. This
  is a true schema signal: a new ``INFO`` field now opens a ``drift:schema`` PR and means
  it. (ClinVar's VCF is sites-only, so ``FORMAT`` is legitimately empty; it is recorded
  anyway so that a future release growing genotype columns is not silently invisible.)
* ``checksums`` -- sha256 over the *structure-declaring* subset of the header only
  (``_SCHEMA_META_PREFIXES``), never the whole header. Hashing the whole header puts
  ``##fileDate`` inside a ``SCHEMA_KEY``, which reproduces #333 exactly: the date moves
  every week, so the tier would be constant every week. This hash is still strictly
  stronger than the ID lists above, because it also catches a ``Number``/``Type``
  changing on a field whose ID did not.
* ``extras.content_length`` -- still compared, so a release still opens a PR, but as
  ``drift:routine`` and batched with the rest.
* ``source_version`` -- ``Last-Modified``. Unlike hgnc, ClinVar does not republish
  byte-identical content under a fresh timestamp, so this identifies a release rather
  than generating noise, and stays in the compared surface.

Fails closed throughout: a partial or unrecognisable header raises ``DriftProbeError``
rather than recording a fingerprint that a later run would compare equal to.
"""

from __future__ import annotations

import hashlib
import io
import zlib
from collections.abc import Mapping
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.clinvar.shared.constants import CLINVAR_FTP_BASE

PROBE_VERSION = 2
_FILENAME = "clinvar.vcf.gz"
_URL = f"{CLINVAR_FTP_BASE}/{_FILENAME}"
_TIMEOUT_S = (5.0, 30.0)

#: Bytes to request. 64 KiB covered the whole header with ~10x margin when measured; the
#: read is hard-capped at this regardless of what the server sends, so a range-ignoring
#: 200 response costs one chunk, not the whole ~194 MB file.
_RANGE_BYTES = 65536


def _fetch_header_bytes() -> tuple[bytes, Mapping[str, str]]:
    """Return the leading bytes of the VCF plus the response headers.

    Streams and stops at ``_RANGE_BYTES``. That cap is the real protection: if NCBI ever
    stops honouring ``Range`` and answers ``200`` with the full body, this reads one
    chunk and closes the connection instead of pulling ~194 MB on a cron.

    Returns ``resp.headers`` itself rather than ``dict(resp.headers)``: requests' own
    ``CaseInsensitiveDict`` is the point. HTTP/2 lower-cases every header name, so a
    plain dict would make ``.get("Content-Range")`` miss and the probe would fail closed
    on a perfectly good response.
    """
    try:
        with requests.get(
            _URL,
            headers={"Range": f"bytes=0-{_RANGE_BYTES - 1}"},
            timeout=_TIMEOUT_S,
            stream=True,
            allow_redirects=True,
        ) as resp:
            # Inside the `with`, so a 4xx/5xx still releases the pooled connection.
            resp.raise_for_status()
            buf = io.BytesIO()
            for chunk in resp.iter_content(8192):
                buf.write(chunk)
                if buf.tell() >= _RANGE_BYTES:
                    break
            return buf.getvalue()[:_RANGE_BYTES], resp.headers
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc


def _decompress_bgzf_prefix(raw: bytes) -> str:
    """Decompress as many whole BGZF members as the prefix contains.

    A truncated final member is expected and fine. Note its partial output is *kept*,
    not discarded -- ``decompress`` runs before the ``eof`` check -- so the returned text
    can end mid-line. Callers must not assume the last line is complete; the ``#CHROM``
    guard in ``fetch_fingerprint`` is what makes that safe, since ``#CHROM`` is terminal
    in a VCF header and every ``##INFO`` line necessarily precedes it intact.

    ``gzip.GzipFile`` is not usable here: it raises ``EOFError`` on the truncation and
    yields nothing at all, rather than the complete members that precede it.
    """
    out = bytearray()
    pos = 0
    while pos < len(raw):
        dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            out += dec.decompress(raw[pos:])
        except zlib.error:
            break
        if not dec.eof:
            break  # truncated trailing member; keep what decompressed
        consumed = len(raw) - pos - len(dec.unused_data)
        if consumed <= 0:
            break
        pos += consumed
    return out.decode("utf-8", errors="replace")


#: Meta lines that describe the file's *structure*. Everything else in the header --
#: ``##fileDate``, ``##source``, and the ``##bcftools_*`` lines a subset carries -- is
#: provenance, and provenance moves on a release that changed nothing structural.
#:
#: Hashing the whole header instead is the trap this probe was written to escape and
#: fell into anyway: ``##fileDate`` changes *every* week, so ``checksums`` moved every
#: week, and since ``checksums`` is a SCHEMA_KEY to the drift bot the tier was constant
#: again -- issue #333 reproduced through a different field. Verified against the live
#: file on 2026-09-21: ``##fileDate=2026-09-13`` sits above ``#CHROM``.
#:
#: ``##reference`` and ``##fileformat`` are structural: a build change or a VCF spec
#: bump both break a builder, and both should open their own PR.
_SCHEMA_META_PREFIXES = (
    "##fileformat=",
    "##reference=",
    "##INFO=",
    "##FORMAT=",
    "##FILTER=",
    "##ALT=",
    "##contig=",
    "#CHROM",
)


def _schema_text(header_text: str) -> str:
    """The structure-declaring subset of the header, in file order.

    Hashing this rather than ``header_text`` is what makes ``checksums`` a schema
    signal. It is strictly stronger than the ``headers`` ID lists on their own: those
    record which ``INFO`` fields exist, while this also catches a ``Number`` or ``Type``
    changing on a field whose ID did not.
    """
    return "\n".join(
        line
        for line in header_text.splitlines()
        if line.startswith(_SCHEMA_META_PREFIXES)
    )


def _parse_meta_ids(header_text: str, kind: str) -> list[str]:
    """IDs declared by ``##<kind>=<ID=...,...>`` lines, sorted and de-duplicated."""
    prefix = f"##{kind}=<ID="
    ids = {
        line[len(prefix) :].split(",", 1)[0].rstrip(">")
        for line in header_text.splitlines()
        if line.startswith(prefix)
    }
    return sorted(ids)


def fetch_fingerprint() -> dict:
    """Fingerprint ClinVar's VCF header (schema) and size (content)."""
    raw, http_headers = _fetch_header_bytes()
    if not raw:
        raise DriftProbeError("ClinVar returned an empty body for the ranged request.")

    text = _decompress_bgzf_prefix(raw)
    lines = text.splitlines()

    # The #CHROM line terminates the header. Without it the prefix was too small or the
    # file is no longer BGZF, and every ID list below would be silently partial -- which
    # would bake an incomplete baseline that later runs compare equal to.
    chrom = [ln for ln in lines if ln.startswith("#CHROM")]
    if not chrom:
        raise DriftProbeError(
            f"no #CHROM line within the first {_RANGE_BYTES} bytes of {_FILENAME} "
            f"({len(raw)} bytes fetched, {len(lines)} header lines decompressed); "
            "the header has grown past the probe's range or the file is not BGZF."
        )

    header_text = text[: text.index(chrom[0]) + len(chrom[0])]
    info_ids = _parse_meta_ids(header_text, "INFO")
    if not info_ids:
        raise DriftProbeError(
            "ClinVar VCF header declared no ##INFO fields; the format has changed."
        )
    format_ids = _parse_meta_ids(header_text, "FORMAT")

    content_length = http_headers.get("Content-Length")
    # On a 206 this is the RANGE length, not the file size -- useless as a content signal
    # and actively misleading, since it would be constant across every release.
    if http_headers.get("Content-Range"):
        content_length = http_headers["Content-Range"].rsplit("/", 1)[-1]
    if not content_length or not content_length.isdigit():
        raise DriftProbeError(
            "could not determine the full size of clinvar.vcf.gz "
            f"(Content-Range={http_headers.get('Content-Range')!r}, "
            f"Content-Length={http_headers.get('Content-Length')!r})"
        )

    return {
        "probe_version": PROBE_VERSION,
        "source_version": http_headers.get("Last-Modified"),
        "headers": {_FILENAME: {"INFO": info_ids, "FORMAT": format_ids}},
        "checksums": {
            _FILENAME: hashlib.sha256(
                _schema_text(header_text).encode("utf-8")
            ).hexdigest()
        },
        "extras": {"content_length": content_length},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
