"""Post-install smoke test for hvantk.

Verifies that Hail initialises correctly, generates a small synthetic
MatrixTable via ``hl.balding_nichols_model``, and prints its schema with
``describe()``.

Usage
-----
    hvantk check-install
"""

import logging
import sys

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


@click.command("check-install", context_settings=CONTEXT_SETTINGS)
def check_install_cmd() -> None:
    """Smoke-test the hvantk installation.

    \b
    Performs the following checks:
      1. Initialises Hail and prints version info.
      2. Generates a small synthetic MatrixTable
         (balding_nichols_model: 3 pops, 50 samples, 200 variants).
      3. Calls describe() on the MatrixTable to print its schema.
    """
    from hvantk.core.hail_context import init_hail

    # -- 1. Hail init ----------------------------------------------------------
    click.echo("1/3  Initialising Hail ...")
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
    except Exception as exc:
        click.secho(f"FAIL  balding_nichols_model failed: {exc}", fg="red")
        sys.exit(1)

    click.echo(f"  OK  {mt.count_rows()} variants x {mt.count_cols()} samples")

    # -- 3. Describe schema ----------------------------------------------------
    click.echo("3/3  MatrixTable schema:")
    mt.describe()

    click.secho("\nAll checks passed.", fg="green")
