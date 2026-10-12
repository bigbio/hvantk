"""GTEx eQTL drift probe: the pinned v11 release, read from the public adult-gtex bucket.

Runs OFFLINE - the Cloud Storage listing and the README are stubbed via requests_mock.
URLs are pinned as literals rather than imported from the probe (issue #341, item 2), so
a probe aimed at the wrong bucket or object raises ``NoMockAddress`` instead of passing.
"""

from __future__ import annotations

import hashlib

import pytest
import requests
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.gtex_eqtl.drift_probe import fetch_fingerprint

LIST_URL = "https://storage.googleapis.com/storage/v1/b/adult-gtex/o"
README_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-qtl/v11/single-tissue-cis-qtl/"
    "README_eQTL_v11.txt"
)
RELEASE_PREFIX = "bulk-qtl/v11/single-tissue-cis-qtl/"

_RULE = "#" + "-" * 78

# The section the builder's input is documented in; the probe must hash exactly this.
_SIGNIF_SECTION = (
    "File extensions: *.signif_pairs.parquet\n"
    "  Column headers:\n"
    "\n"
    "    phenotype_id:             Phenotype ID\n"
    "    variant_id:               variant ID\n"
    "    af:                       allele frequency of the ALT allele (in-sample)"
)

_README = "\n".join(
    [
        "#" + "=" * 78,
        "#  README for cis-QTLs for the v11 GTEx release",
        "#" + "=" * 78,
        "",
        _RULE,
        "#  eGenes / sGenes / apaGenes",
        _RULE,
        "",
        "  File extension: *.eGenes.txt.gz",
        "  Column headers:",
        "",
        "    gene_name:                GENCODE gene name",
        "",
        _RULE,
        "#  Significant variant-gene pairs",
        _RULE,
        "",
        "  " + _SIGNIF_SECTION,
        "",
        "",
        _RULE,
        "#  All variant-gene pairs tested",
        _RULE,
        "",
        "  File extension: *.allpairs.txt.gz",
    ]
)

_ARCHIVE = {"md5Hash": "O68VeoUuGIDxrhyAzSwQHw==", "size": "3288535040"}
_README_OBJECT = {"md5Hash": "tT4pzpwctWjgqXHqiTd2Aw==", "size": "7182"}


def _mock_bucket(m, *, objects=("GTEx_Analysis_v11_eQTL.tar", "README_eQTL_v11.txt")):
    """Answer both listings the probe makes, keyed on the requested prefix."""

    def listing(request, context):
        prefix = request.qs["prefix"][0]
        if prefix == "bulk-qtl/":
            return {"prefixes": ["bulk-qtl/v8/", "bulk-qtl/v11/", "bulk-qtl/v10/"]}
        assert prefix == RELEASE_PREFIX, prefix
        return {
            "items": [
                {
                    "name": RELEASE_PREFIX + name,
                    **(_ARCHIVE if name.endswith(".tar") else _README_OBJECT),
                }
                for name in objects
            ]
        }

    m.get(LIST_URL, json=listing)


def test_fetch_fingerprint_shape():
    with requests_mock.Mocker() as m:
        _mock_bucket(m)
        m.get(README_URL, text=_README)
        fp = fetch_fingerprint()

    readme_key = "README_eQTL_v11.txt:signif_pairs"
    assert fp["probe_version"] == 2
    assert fp["source_version"] == "v11"
    assert fp["headers"] == {
        "gs://adult-gtex/bulk-qtl/": ["v10", "v11", "v8"],
        "gs://adult-gtex/bulk-qtl/v11/single-tissue-cis-qtl/": [
            "GTEx_Analysis_v11_eQTL.tar",
            "README_eQTL_v11.txt",
        ],
        # Only the signif_pairs section: `gene_name` lives in the eGenes section.
        readme_key: ["phenotype_id", "variant_id", "af"],
    }
    assert fp["checksums"] == {
        readme_key: hashlib.sha256(_SIGNIF_SECTION.encode("utf-8")).hexdigest()
    }
    assert fp["extras"] == {"GTEx_Analysis_v11_eQTL.tar": _ARCHIVE}
    assert "fetched_at" in fp


def test_a_missing_archive_is_layout_drift_not_a_failure():
    """The folder listing under `headers` reports it; extras records the archive as null."""
    with requests_mock.Mocker() as m:
        _mock_bucket(m, objects=("README_eQTL_v11.txt",))
        m.get(README_URL, text=_README)
        fp = fetch_fingerprint()
    assert fp["headers"]["gs://adult-gtex/bulk-qtl/v11/single-tissue-cis-qtl/"] == [
        "README_eQTL_v11.txt"
    ]
    assert fp["extras"] == {"GTEx_Analysis_v11_eQTL.tar": None}


def test_fetch_fingerprint_raises_when_readme_drops_signif_pairs():
    readme = _README.replace("*.signif_pairs.parquet", "*.pairs.txt.gz")
    with requests_mock.Mocker() as m:
        _mock_bucket(m)
        m.get(README_URL, text=readme)
        with pytest.raises(DriftProbeError, match="signif_pairs"):
            fetch_fingerprint()


def test_fetch_fingerprint_raises_on_request_failure():
    with requests_mock.Mocker() as m:
        m.get(LIST_URL, exc=requests.ConnectionError("boom"))
        with pytest.raises(DriftProbeError):
            fetch_fingerprint()
