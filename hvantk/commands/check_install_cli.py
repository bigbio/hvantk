"""Post-install smoke test for hvantk.

Verifies that Hail initialises correctly, generates a small synthetic
MatrixTable via ``hl.balding_nichols_model``, and prints its schema with
``describe()``.

Usage
-----
    hvantk check-install
"""

import logging
import os
import sys

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)

_PROXY_VARS = ("HTTP_PROXY", "HTTPS_PROXY", "http_proxy", "https_proxy")


def _check_proxy() -> bool:
    """Warn if proxy env vars are set and NO_PROXY does not cover localhost.

    Returns True if a problematic proxy was detected.
    """
    active = {k: v for k in _PROXY_VARS if (v := os.environ.get(k))}
    if not active:
        return False

    no_proxy = os.environ.get("NO_PROXY", "") + "," + os.environ.get("no_proxy", "")
    localhost_covered = any(
        tok in no_proxy for tok in ("localhost", "127.0.0.1", "0.0.0.0")
    )
    if localhost_covered:
        return False

    click.secho("WARN  HTTP proxy detected without localhost bypass:", fg="yellow")
    for k, v in active.items():
        click.secho(f"        {k}={v}", fg="yellow")
    click.secho(
        "      Hail communicates with Spark over localhost HTTP.\n"
        "      A proxy intercepting this traffic causes orjson.JSONDecodeError.\n"
        "      Fix: unset HTTP_PROXY HTTPS_PROXY http_proxy https_proxy\n"
        "           export NO_PROXY='localhost,127.0.0.1,0.0.0.0,::1'\n"
        "      See: https://github.com/bigbio/hvantk/issues/36",
        fg="yellow",
    )
    return True


@click.command("check-install", context_settings=CONTEXT_SETTINGS)
def check_install_cmd() -> None:
    """Smoke-test the hvantk installation.

    \b
    Performs the following checks:
      0. Checks for HTTP proxy that may break Hail (see issue #36).
      1. Initialises Hail and prints version info.
      2. Generates a small synthetic MatrixTable
         (balding_nichols_model: 3 pops, 50 samples, 200 variants).
      3. Calls describe() on the MatrixTable to print its schema.
    """
    # -- 0. Proxy check --------------------------------------------------------
    click.echo("0/3  Checking environment ...")
    proxy_found = _check_proxy()
    if proxy_found:
        click.secho(
            "      Continuing anyway — Hail operations will likely fail.\n",
            fg="yellow",
        )
    else:
        click.echo("  OK  No problematic proxy detected")

    # -- 1. Hail init ----------------------------------------------------------
    click.echo("1/3  Initialising Hail ...")
    from hvantk.core.hail_context import init_hail

    try:
        init_hail()
    except Exception as exc:
        click.secho(f"FAIL  Hail initialisation failed: {exc}", fg="red")
        sys.exit(1)

    import hail as hl

    click.echo(f"  OK  Hail version {hl.version()}")

    # -- 2. Generate synthetic MatrixTable -------------------------------------
    click.echo("2/3  Generating synthetic MatrixTable (balding_nichols_model) ...")
    try:
        mt = hl.balding_nichols_model(
            n_populations=3,
            n_samples=50,
            n_variants=200,
            n_partitions=2,
        )
        n_rows, n_cols = mt.count()
    except Exception as exc:
        if "orjson" in str(type(exc).__name__).lower() or "orjson" in str(exc).lower():
            click.secho(f"FAIL  {exc}", fg="red")
            click.secho(
                "\n      This is likely caused by an HTTP proxy intercepting "
                "Hail's localhost traffic.\n"
                "      Run:  unset HTTP_PROXY HTTPS_PROXY http_proxy https_proxy\n"
                "            export NO_PROXY='localhost,127.0.0.1,0.0.0.0,::1'\n"
                "      See:  https://github.com/bigbio/hvantk/issues/36",
                fg="yellow",
            )
        else:
            click.secho(f"FAIL  balding_nichols_model failed: {exc}", fg="red")
        sys.exit(1)

    click.echo(f"  OK  {n_rows} variants x {n_cols} samples")

    # -- 3. Describe schema ----------------------------------------------------
    click.echo("3/3  MatrixTable schema:")
    mt.describe()

    click.secho("\nAll checks passed.", fg="green")
