"""ClinVar drift probe: schema signal from the VCF header, content signal from size.

Runs OFFLINE - the ranged GET is stubbed with requests_mock over a synthetic BGZF
payload, so CI never touches NCBI.

These tests exist because ``probe_version`` 1 could not fail usefully: it recorded
``Content-Length`` under ``headers``, a *schema* key, so every ClinVar release tiered
``drift:schema`` and the tier stopped carrying information (issue #333). The important
assertions here are therefore not "the shape is right" but "the schema signal and the
content signal are in different places, and each moves only when it should".
"""

from __future__ import annotations

import zlib

import pytest
import requests
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.clinvar.drift_probe import fetch_fingerprint, _RANGE_BYTES, _URL

_HEADER = (
    "##fileformat=VCFv4.1\n"
    '##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="the ClinVar Allele ID">\n'
    '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">\n'
    '##INFO=<ID=AF_ESP,Number=1,Type=Float,Description="allele frequencies from GO-ESP">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "1\t100\t.\tA\tG\t.\t.\tALLELEID=1\n"
)


def _bgzf(payload: str, members: int = 3) -> bytes:
    """Concatenated gzip members, as BGZF is -- not a single deflate stream.

    Splitting across members is the point: the probe must reassemble several of them,
    which is what a real ranged read returns.
    """
    chunks, out = payload.encode(), b""
    step = max(1, len(chunks) // members)
    for i in range(0, len(chunks), step):
        co = zlib.compressobj(9, zlib.DEFLATED, 16 + zlib.MAX_WBITS)
        out += co.compress(chunks[i : i + step]) + co.flush()
    return out


def _stub(m, body: bytes, *, full_size: str = "193427450", **extra):
    headers = {
        "Last-Modified": "Wed, 01 Jan 2026 00:00:00 GMT",
        "Content-Range": f"bytes 0-{len(body) - 1}/{full_size}",
        "Content-Length": str(len(body)),
    }
    headers.update(extra)
    m.get(_URL, content=body, headers=headers, status_code=206)


def test_schema_signal_is_the_real_info_list():
    """``headers`` must carry VCF field IDs -- not HTTP metadata, as v1 did."""
    with requests_mock.Mocker() as m:
        _stub(m, _bgzf(_HEADER))
        fp = fetch_fingerprint()

    assert fp["probe_version"] == 2
    assert fp["source_version"] == "Wed, 01 Jan 2026 00:00:00 GMT"
    assert fp["headers"]["clinvar.vcf.gz"]["INFO"] == ["AF_ESP", "ALLELEID", "CLNSIG"]
    # Sites-only today, but recorded so a future genotype column is not invisible.
    assert fp["headers"]["clinvar.vcf.gz"]["FORMAT"] == []
    assert fp["checksums"]["clinvar.vcf.gz"]
    assert "fetched_at" in fp


def test_content_length_is_the_full_file_not_the_range():
    """A 206 reports the *range* length; recording that would freeze the signal.

    ``Content-Length`` on a partial response describes the slice, which is constant
    across releases -- so the content signal has to come from ``Content-Range``.
    """
    with requests_mock.Mocker() as m:
        _stub(m, _bgzf(_HEADER), full_size="193427450")
        fp = fetch_fingerprint()
    assert fp["extras"]["content_length"] == "193427450"


def test_size_change_alone_leaves_the_schema_signal_untouched():
    """ClinVar adding variants must not move ``headers``/``checksums`` (#333).

    This is the whole point of the rewrite: the drift bot tiers on those two keys, so if
    a size change perturbs either, every release is a schema change again.
    """
    with requests_mock.Mocker() as m:
        _stub(m, _bgzf(_HEADER), full_size="193427450")
        before = fetch_fingerprint()
    with requests_mock.Mocker() as m:
        _stub(m, _bgzf(_HEADER), full_size="199999999")
        after = fetch_fingerprint()

    assert before["headers"] == after["headers"]
    assert before["checksums"] == after["checksums"]
    assert before["extras"]["content_length"] != after["extras"]["content_length"]


def test_a_new_info_field_moves_both_schema_keys():
    with requests_mock.Mocker() as m:
        _stub(m, _bgzf(_HEADER))
        before = fetch_fingerprint()
    grown = _HEADER.replace(
        "#CHROM",
        '##INFO=<ID=ONCDN,Number=.,Type=String,Description="new">\n#CHROM',
    )
    with requests_mock.Mocker() as m:
        _stub(m, _bgzf(grown))
        after = fetch_fingerprint()

    assert "ONCDN" in after["headers"]["clinvar.vcf.gz"]["INFO"]
    assert before["checksums"] != after["checksums"]


def test_header_beyond_the_fetched_range_fails_closed():
    """A truncated header must raise, never record a partial ID list.

    Recording it would bake an incomplete baseline that every later run compares equal
    to -- the probe would report clean having quietly stopped seeing the schema.
    """
    headless = (
        "##fileformat=VCFv4.1\n"
        + '##INFO=<ID=X,Number=1,Type=Integer,Description="x">\n'
    )
    with requests_mock.Mocker() as m:
        _stub(m, _bgzf(headless))
        with pytest.raises(DriftProbeError, match="#CHROM"):
            fetch_fingerprint()


def test_header_without_info_fields_fails_closed():
    with requests_mock.Mocker() as m:
        _stub(m, _bgzf("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\n"))
        with pytest.raises(DriftProbeError, match="no ##INFO"):
            fetch_fingerprint()


def test_unknown_full_size_fails_closed():
    """Without ``Content-Range`` the content signal cannot be trusted, so refuse."""
    with requests_mock.Mocker() as m:
        body = _bgzf(_HEADER)
        m.get(
            _URL,
            content=body,
            headers={"Last-Modified": "Wed, 01 Jan 2026 00:00:00 GMT"},
            status_code=200,
        )
        with pytest.raises(DriftProbeError, match="could not determine the full size"):
            fetch_fingerprint()


def test_request_failure_raises():
    with requests_mock.Mocker() as m:
        m.get(_URL, exc=requests.ConnectionError("boom"))
        with pytest.raises(DriftProbeError, match="HTTP failure"):
            fetch_fingerprint()


def test_range_header_is_actually_sent():
    """The cap is what stops a range-ignoring server costing a ~500 MB download."""
    with requests_mock.Mocker() as m:
        _stub(m, _bgzf(_HEADER))
        fetch_fingerprint()
        assert m.last_request.headers["Range"] == f"bytes=0-{_RANGE_BYTES - 1}"
