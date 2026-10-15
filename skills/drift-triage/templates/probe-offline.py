"""Answer hvantk's drift probes from a local file instead of the network.

Coverage: probes that call the module-level ``requests.head`` / ``requests.get``
functions directly -- ``hvantk/skills/clingen/drift_probe.py`` and
``hvantk/skills/hgnc/drift_probe.py`` both do. A probe built on ``requests.Session`` or
any other HTTP client is NOT covered: this module patches only the two module-level
functions, and a probe that never calls them observes nothing different.

Needs ``hvantk`` importable (run inside the container image, or with the checkout on
``PYTHONPATH``); no Hail session is started -- the clingen/hgnc probes this covers need
none, and nothing here calls ``init_hail()``.

Three subcommands:

    probe-offline.py fingerprint --dataset PROVIDER:DATASET --file PATH
        Print the probe's fingerprint for PATH, as if PATH were served at the probe's
        own URL.

    probe-offline.py diff --dataset PROVIDER:DATASET --file NEW --against OLD
        Print, per header filename key, the columns NEW adds and removes relative to
        OLD -- the {added, removed} shape of a drift diff
        (hvantk.core.plugin.drift_runner._compare_fingerprints), restricted to
        `headers`. NEW is "observed", OLD is "expected", matching that same convention.

    probe-offline.py drift --dataset PROVIDER:DATASET --file PATH -- [hvantk drift options]
        Run `hvantk drift PROVIDER:DATASET [options]` in-process with upstream answered
        from PATH, and exit with the CLI's own exit code (hvantk/tools/plugins/drift_cli.py:
        0 clean, 1 drifted, 2 probe failed, 3 registry error, 4 usage error). The
        pass-through options must follow a literal `--`.

Used where the real probe cannot reach upstream (a sandboxed eval, a login node or
compute node with no egress) or where the signal must be pinned to one specific file
rather than whatever upstream currently serves.
"""

from __future__ import annotations

import argparse
import contextlib
import email.utils
import json
import sys
from pathlib import Path


@contextlib.contextmanager
def _answer_from_file(path):
    """Patch ``requests.head``/``requests.get`` so every call answers from ``path``.

    Builds a fresh ``requests.Response`` per call: status 200, headers carrying
    Content-Length / Last-Modified / Content-Type for ``path``, and body bytes read
    from ``path`` (empty for HEAD, matching a real HEAD response). Restores both
    functions on exit, even if the probe raises.
    """
    import requests
    from requests.structures import CaseInsensitiveDict

    file_path = Path(path)
    data = file_path.read_bytes()
    last_modified = email.utils.formatdate(file_path.stat().st_mtime, usegmt=True)

    def _make_response(url, body):
        resp = requests.Response()
        resp.status_code = 200
        resp.reason = "OK"
        resp.url = url
        resp.headers = CaseInsensitiveDict(
            {
                "Content-Length": str(len(data)),
                "Last-Modified": last_modified,
                "Content-Type": "text/plain",
            }
        )
        resp._content = body
        resp._content_consumed = True
        resp.encoding = "utf-8"
        resp.raw = None
        return resp

    def fake_head(url, *args, **kwargs):
        return _make_response(url, b"")

    def fake_get(url, *args, **kwargs):
        return _make_response(url, data)

    orig_head, orig_get = requests.head, requests.get
    requests.head, requests.get = fake_head, fake_get
    try:
        yield
    finally:
        requests.head, requests.get = orig_head, orig_get


def _dataset_probe(dataset):
    """The dataset's drift-probe callable, resolved through the real plugin registry."""
    from hvantk.core.plugin.loader import get_registry

    return get_registry().get_dataset(dataset).drift_probe


def cmd_fingerprint(args) -> int:
    probe = _dataset_probe(args.dataset)
    with _answer_from_file(args.file):
        fingerprint = probe()
    print(json.dumps(fingerprint, indent=2, default=str))
    return 0


def cmd_diff(args) -> int:
    probe = _dataset_probe(args.dataset)
    with _answer_from_file(args.file):
        new_fp = probe()
    with _answer_from_file(args.against):
        old_fp = probe()
    new_headers = new_fp.get("headers") or {}
    old_headers = old_fp.get("headers") or {}
    out = {}
    for name in sorted(set(new_headers) | set(old_headers)):
        new_cols = new_headers.get(name, [])
        old_cols = old_headers.get(name, [])
        out[name] = {
            "added": [c for c in new_cols if c not in old_cols],
            "removed": [c for c in old_cols if c not in new_cols],
        }
    print(json.dumps(out, indent=2))
    return 0


def cmd_drift(args) -> int:
    from hvantk.hvantk import cli

    drift_args = list(args.drift_args)
    if drift_args and drift_args[0] == "--":
        drift_args = drift_args[1:]

    with _answer_from_file(args.file):
        try:
            cli.main(
                args=["drift", args.dataset, *drift_args],
                prog_name="hvantk",
                standalone_mode=True,
            )
        except SystemExit as exc:
            code = exc.code
            if isinstance(code, int):
                return code
            return 0 if code is None else 1
    # cli.main(standalone_mode=True) always ends by calling sys.exit (see
    # hvantk/tools/plugins/drift_cli.py), so the SystemExit above is the only return
    # path in practice; this is defensive, not expected to run.
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="command", required=True)

    fp = sub.add_parser("fingerprint", help="print a probe's fingerprint for a file")
    fp.add_argument("--dataset", required=True, metavar="PROVIDER:DATASET")
    fp.add_argument("--file", required=True)
    fp.set_defaults(fn=cmd_fingerprint)

    df = sub.add_parser("diff", help="diff two files' header columns through a probe")
    df.add_argument("--dataset", required=True, metavar="PROVIDER:DATASET")
    df.add_argument("--file", required=True, help='the new file ("observed")')
    df.add_argument("--against", required=True, help='the old file ("expected")')
    df.set_defaults(fn=cmd_diff)

    dr = sub.add_parser(
        "drift", help="run `hvantk drift` with upstream answered from a file"
    )
    dr.add_argument("--dataset", required=True, metavar="PROVIDER:DATASET")
    dr.add_argument("--file", required=True)
    dr.add_argument(
        "drift_args",
        nargs=argparse.REMAINDER,
        metavar="-- [hvantk drift options]",
        help="pass-through options for `hvantk drift`, after a literal --",
    )
    dr.set_defaults(fn=cmd_drift)

    args = ap.parse_args()
    try:
        return args.fn(args)
    except Exception as exc:  # noqa: BLE001 - a CLI entry point, not a library call
        print(f"error: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
