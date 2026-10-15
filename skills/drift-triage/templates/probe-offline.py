"""Answer hvantk's drift probes from a local file instead of the network.

Coverage: probes that call the module-level ``requests.head``, ``requests.get`` or
``requests.request`` functions. clingen and hgnc call ``head``/``get`` directly; gencc
goes through ``hvantk.core.utils.http.request_with_retry``, which -- called with no
``session=`` argument, as the gencc probe does -- calls module-level
``requests.request``, so it is served too. A probe built on its OWN
``requests.Session`` (gnomad-metrics) is NOT served: for the duration of the context
below, ``requests.Session.send`` is patched to raise
``RuntimeError("probe-offline: the probe attempted a network call the offline form "
"does not serve")`` instead, so such a probe fails loudly -- as a ``DriftProbeError``,
which ``hvantk drift`` reports as ``probe_failed`` -- rather than silently reaching the
real network or answering for the wrong file.

Needs ``hvantk`` importable (run inside the container image, or with the checkout on
``PYTHONPATH``); no Hail session is started -- the probes this covers need none, and
nothing here calls ``init_hail()``.

Three subcommands:

    probe-offline.py fingerprint --dataset PROVIDER:DATASET --file PATH
        Print the probe's fingerprint for PATH, as if PATH were served at the probe's
        own URL. Exits 5 if the probe made no head/get/request call at all (Session-based
        or otherwise uncovered) -- there is nothing to print.

    probe-offline.py diff --dataset PROVIDER:DATASET --file NEW (--against OLD | --against-committed)
        Print, per header filename key, the COLUMNS NEW adds and removes relative to
        OLD. This is narrower than a full fingerprint's own top-level ``added``/
        ``removed`` (hvantk.core.plugin.drift_runner._compare_fingerprints, which the
        `drift` subcommand below surfaces unabridged) -- here it is always just the
        `headers` entry's column list. NEW is "observed", OLD is "expected", matching
        that same convention. Exactly one of ``--against``/``--against-committed`` is
        required: ``--against FILE`` probes a second local file the same way as NEW;
        ``--against-committed`` skips probing anything for OLD and instead reads the
        dataset's own committed ``drift_fingerprint.json`` (via
        ``get_registry().get_dataset(ds).test_paths.drift_fingerprint``) and uses ITS
        ``headers`` entry -- the real baseline a repair must match, as opposed to a
        plugin's test fixture, which can carry far fewer columns than upstream ever
        did. Exits 5 under the same condition as `fingerprint` (for NEW; the committed
        fingerprint is read directly, not probed, so it has no such failure mode).

    probe-offline.py drift --dataset PROVIDER:DATASET --file PATH -- [hvantk drift options]
        Run `hvantk drift PROVIDER:DATASET [options]` in-process with upstream answered
        from PATH, and exit with the CLI's own exit code (hvantk/tools/plugins/drift_cli.py:
        0 clean, 1 drifted, 2 probe failed -- including a Session-blocked probe --,
        3 registry error, 4 usage error). The pass-through options must follow a
        literal `--`. An argparse error on THIS template's own arguments (a missing
        --dataset, an unknown subcommand) exits 2, like any other argparse usage error.

Any other failure in this template itself (an unreadable file, an unknown dataset,
anything not covered by the two paths above) exits 5.

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


class _ServedCalls:
    """How many head/get/request calls `_answer_from_file` actually served."""

    def __init__(self):
        self.count = 0


_BLOCKED_MESSAGE = (
    "probe-offline: the probe attempted a network call the offline form does not serve"
)


@contextlib.contextmanager
def _answer_from_file(path):
    """Patch ``requests.head``/``get``/``request`` so every call answers from ``path``,
    and block anything going through a ``requests.Session`` instead.

    Builds a fresh ``requests.Response`` per call: status 200, headers carrying
    Content-Length / Last-Modified / Content-Type for ``path``, and body bytes read
    from ``path`` (empty for HEAD, matching a real HEAD response). Yields a
    `_ServedCalls` counter so the caller can tell whether the probe actually used one
    of the three served entry points. Restores everything on exit, even if the probe
    raises.
    """
    import requests
    from requests.structures import CaseInsensitiveDict

    file_path = Path(path)
    data = file_path.read_bytes()
    last_modified = email.utils.formatdate(file_path.stat().st_mtime, usegmt=True)
    served = _ServedCalls()

    def _make_response(url, body):
        served.count += 1
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

    def fake_request(method, url, *args, **kwargs):
        body = b"" if str(method).upper() == "HEAD" else data
        return _make_response(url, body)

    def _blocked_send(self, request, **kwargs):
        raise RuntimeError(_BLOCKED_MESSAGE)

    orig_head, orig_get, orig_request = requests.head, requests.get, requests.request
    orig_send = requests.Session.send
    requests.head, requests.get, requests.request = fake_head, fake_get, fake_request
    requests.Session.send = _blocked_send
    try:
        yield served
    finally:
        requests.head, requests.get, requests.request = (
            orig_head,
            orig_get,
            orig_request,
        )
        requests.Session.send = orig_send


def _dataset_probe(dataset):
    """The dataset's drift-probe callable, resolved through the real plugin registry."""
    from hvantk.core.plugin.loader import get_registry

    return get_registry().get_dataset(dataset).drift_probe


def _not_served_error(dataset: str) -> None:
    print(
        f"error: the offline form served no requests.head/get/request call for "
        f"{dataset}'s probe; it may use requests.Session (not covered here) or "
        "another HTTP client entirely",
        file=sys.stderr,
    )


def cmd_fingerprint(args) -> int:
    probe = _dataset_probe(args.dataset)
    with _answer_from_file(args.file) as served:
        fingerprint = probe()
    if served.count == 0:
        _not_served_error(args.dataset)
        return 5
    print(json.dumps(fingerprint, indent=2, default=str))
    return 0


def cmd_diff(args) -> int:
    probe = _dataset_probe(args.dataset)
    with _answer_from_file(args.file) as new_served:
        new_fp = probe()
    if new_served.count == 0:
        _not_served_error(args.dataset)
        return 5
    new_headers = new_fp.get("headers") or {}

    if args.against_committed:
        from hvantk.core.plugin.loader import get_registry

        fp_path = Path(
            get_registry().get_dataset(args.dataset).test_paths.drift_fingerprint
        )
        old_fp = json.loads(fp_path.read_text())
    else:
        with _answer_from_file(args.against) as old_served:
            old_fp = probe()
        if old_served.count == 0:
            _not_served_error(args.dataset)
            return 5
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
    df_against = df.add_mutually_exclusive_group(required=True)
    df_against.add_argument("--against", help='the old file ("expected")')
    df_against.add_argument(
        "--against-committed",
        action="store_true",
        help="compare against the dataset's committed drift_fingerprint.json headers",
    )
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
        return 5


if __name__ == "__main__":
    sys.exit(main())
