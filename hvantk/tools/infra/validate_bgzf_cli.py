import os
import sys
import time
from glob import glob

import click

from hvantk.core.config import CONTEXT_SETTINGS
from hvantk.data.file_utils import validate_bgzf


@click.command("validate-bgzf", context_settings=CONTEXT_SETTINGS)
@click.argument("path", type=click.Path(exists=True))
@click.option(
    "--quick", is_flag=True, help="Header-only check, skip CRC/ISIZE verification"
)
def validate_bgzf_cmd(path, quick):
    """Validate BGZF block integrity of BGZF-compressed files."""
    full = not quick
    mode = "header-only" if quick else "full (header + CRC32 + ISIZE)"

    if os.path.isdir(path):
        vcf_files = sorted(glob(os.path.join(path, "*.vcf.gz")))
        if not vcf_files:
            click.echo(f"No .vcf.gz files found in {path}")
            sys.exit(1)
    else:
        vcf_files = [path]

    click.echo(f"Checking {len(vcf_files)} file(s)")
    click.echo(f"Mode: {mode}")
    click.echo("=" * 80)

    passed = 0
    failed = 0
    for vcf in vcf_files:
        basename = os.path.basename(vcf)
        size_gb = os.path.getsize(vcf) / (1024**3)
        click.echo(f"  {basename} ({size_gb:.1f} GB) ... ", nl=False)

        t0 = time.time()
        is_valid, num_blocks, message = validate_bgzf(vcf, full=full)
        elapsed = time.time() - t0

        if is_valid:
            click.echo(f"OK  ({num_blocks:,} blocks, {elapsed:.1f}s)")
            passed += 1
        else:
            click.echo(f"FAIL at block {num_blocks:,}")
            click.echo(f"    {message}")
            failed += 1

    click.echo("=" * 80)
    click.echo(f"Results: {passed} passed, {failed} failed, {len(vcf_files)} total")

    if failed > 0:
        sys.exit(1)
