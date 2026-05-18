"""`hvantk reprocess <provider:dataset>` orchestrates download -> parse -> build.

The command resolves a dataset from the plugin registry and invokes the
lifecycle stages declared in its plugin.yaml. Operators can resume mid-pipeline
with the `--skip-*` flags.
"""

from __future__ import annotations

import json

import click


@click.command(name="reprocess")
@click.argument("dataset")
@click.option(
    "--raw-dir",
    required=True,
    type=click.Path(file_okay=False),
    help="Where to write raw downloaded files",
)
@click.option(
    "--intermediate",
    required=False,
    type=click.Path(),
    default=None,
    help="Path for parsed intermediate (required when the plugin declares lifecycle.parse)",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(),
    help="Final builder output path",
)
@click.option("--skip-download", is_flag=True, help="Skip the download stage")
@click.option("--skip-parse", is_flag=True, help="Skip the parse stage")
@click.option("--skip-build", is_flag=True, help="Skip the build stage")
@click.option(
    "--check-drift",
    is_flag=True,
    default=True,
    help="Run drift probe after build (default: on)",
)
@click.option(
    "--no-check-drift",
    is_flag=True,
    help="Disable the post-build drift probe",
)
def reprocess_cmd(
    dataset,
    raw_dir,
    intermediate,
    output,
    skip_download,
    skip_parse,
    skip_build,
    check_drift,
    no_check_drift,
):
    """Run download -> parse -> build for a plugin dataset."""
    from hvantk.core import drift_runner, plugin_loader

    reg = plugin_loader.get_registry()
    try:
        spec = reg.get_dataset(dataset)
    except KeyError:
        raise click.ClickException(f"unknown dataset: {dataset}")

    # 1. Download stage
    if not skip_download:
        if spec.download_fn is None:
            raise click.UsageError(
                f"{dataset} has no lifecycle.download declared; "
                "use --skip-download or add it to plugin.yaml"
            )
        click.echo(f"download: {dataset} -> {raw_dir}")
        spec.download_fn(raw_dir=raw_dir)

    # 2. Parse stage
    parsed_path = intermediate
    if not skip_parse and spec.parse_fn is not None:
        if intermediate is None:
            raise click.UsageError(
                "--intermediate is required when the plugin declares lifecycle.parse"
            )
        click.echo(f"parse: {raw_dir} -> {intermediate}")
        spec.parse_fn(raw_dir=raw_dir, output_path=intermediate)
        parsed_path = intermediate
    elif not skip_parse:
        # No parse stage declared: the builder consumes the raw dir directly.
        parsed_path = raw_dir

    # 3. Build stage
    if not skip_build:
        click.echo(f"build: {parsed_path} -> {output}")
        spec.builder(parsed_path, output)

    # 4. Optional drift check
    if check_drift and not no_check_drift:
        result = drift_runner.run_drift_check(dataset)
        click.echo(f"drift: {result.status}")
        if result.status == "drifted" and result.diff:
            click.echo(json.dumps(result.diff, indent=2, default=str))
