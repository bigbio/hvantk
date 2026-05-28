"""`hvantk reprocess <provider:dataset>` orchestrates download -> parse -> build.

The command resolves a dataset from the plugin registry and invokes the
lifecycle stages declared in its plugin.yaml. Operators can resume mid-pipeline
with the `--skip-*` flags.
"""

from __future__ import annotations

import json
import re
from typing import Any

import click


_INT_RE = re.compile(r"^-?\d+$")
_FLOAT_RE = re.compile(
    r"^-?(\d+\.\d*([eE][+-]?\d+)?|\.\d+([eE][+-]?\d+)?|\d+[eE][+-]?\d+)$"
)


def _coerce_plugin_arg_value(value: str) -> Any:
    """Coerce a --plugin-arg string to bool / int / float when it matches.

    Order: bool (``true``/``false``, case-insensitive) -> int -> float -> str.
    Without coercion, builders that declare typed kwargs receive surprising
    values — e.g. ``overwrite=false`` arrives as the truthy string ``"false"``,
    and ``p_threshold=5e-8`` arrives as the string ``"5e-8"`` which breaks
    numeric comparisons.

    Pattern matching is used (rather than ``try/except ValueError``) so that
    the function has a single exit path per branch and does not swallow
    exceptions silently.
    """
    # List coercion: if the raw value contains a comma, treat it as a
    # comma-separated list and recursively coerce each element. Builders
    # that declare list-typed params (e.g. chromosomes: list[str]) thus
    # receive a real list instead of the literal "chr1,chr2,chrX" string.
    # Trade-off: a single element containing a comma cannot be expressed
    # via --plugin-arg under this scheme; in practice plugin args don't
    # carry commas in single values (issue #119, Option A).
    if "," in value:
        return [
            _coerce_plugin_arg_value(part.strip())
            for part in value.split(",")
            if part.strip()
        ]
    low = value.lower()
    if low == "true":
        return True
    if low == "false":
        return False
    if _INT_RE.match(value):
        return int(value)
    if _FLOAT_RE.match(value):
        return float(value)
    return value


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
    "--plugin-arg",
    "plugin_args",
    multiple=True,
    metavar="KEY=VALUE",
    help=(
        "Plugin-specific kwarg forwarded to the lifecycle download, parse, "
        "and Phase B build stages (e.g. --plugin-arg cancer_type=brca for "
        "cptac:phospho, --plugin-arg reference_genome=GRCh38 for clinvar). "
        "Repeatable. Values are coerced: 'true'/'false' -> bool, integer "
        "strings -> int, decimal/scientific -> float, otherwise str."
    ),
)
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
    plugin_args,
    check_drift,
    no_check_drift,
):
    """Run download -> parse -> build for a plugin dataset."""
    from hvantk.core.plugin import drift_runner, loader as plugin_loader

    reg = plugin_loader.get_registry()
    try:
        spec = reg.get_dataset(dataset)
    except KeyError:
        raise click.ClickException(f"unknown dataset: {dataset}")

    # Parse --plugin-arg KEY=VALUE pairs into a kwargs dict forwarded to
    # download_fn, parse_fn, and the Phase B builder. Lets operators target
    # plugin-specific options (e.g. cancer_type for cptac:phospho,
    # reference_genome for clinvar, p_threshold for gtex-eqtl) without
    # per-plugin reprocess flags.
    extras: dict[str, Any] = {}
    for kv in plugin_args:
        if "=" not in kv:
            raise click.UsageError(
                f"--plugin-arg must be KEY=VALUE; got: {kv!r}"
            )
        key, value = kv.split("=", 1)
        if not key:
            raise click.UsageError(
                f"--plugin-arg key must be non-empty; got: {kv!r}"
            )
        extras[key] = _coerce_plugin_arg_value(value)

    # 1. Download stage
    if not skip_download:
        if spec.download_fn is None:
            raise click.UsageError(
                f"{dataset} has no lifecycle.download declared; "
                "use --skip-download or add it to plugin.yaml"
            )
        click.echo(f"download: {dataset} -> {raw_dir}")
        spec.download_fn(raw_dir=raw_dir, **extras)

    # 2. Parse stage
    if not skip_parse and spec.parse_fn is not None:
        if intermediate is None:
            raise click.UsageError(
                "--intermediate is required when the plugin declares lifecycle.parse"
            )
        click.echo(f"parse: {raw_dir} -> {intermediate}")
        spec.parse_fn(raw_dir=raw_dir, output_path=intermediate, **extras)
        parsed_path = intermediate
    elif not skip_parse:
        # No parse stage declared: the builder consumes the raw dir directly.
        parsed_path = raw_dir
    else:
        # skip_parse=True: prefer --intermediate if given, else the raw_dir.
        # The prior version left parsed_path as None here, which then crashed
        # the build stage with a confusing TypeError deep inside the builder.
        parsed_path = intermediate if intermediate is not None else raw_dir

    # 3. Build stage
    if not skip_build:
        click.echo(f"build: {parsed_path} -> {output}")
        if spec.artifact_type is None:
            # Legacy (Phase A) plugin not yet migrated to Phase B contract.
            # Fall back to the old shape: spec.builder(input, output) writes
            # directly to disk. No provenance / artifact_type validation.
            spec.builder(parsed_path, output)
        else:
            # Phase B contract: orchestrator handles BuildContext, validation,
            # and save. Returns the stamped Provenance.
            from hvantk.core.plugin.run_builder import run_builder_for_spec
            from pathlib import Path
            run_builder_for_spec(
                spec,
                parsed_input=parsed_path,
                output_path=Path(output),
                plugin_version=spec.plugin_version or "<unknown>",
                **extras,
            )

    # 4. Optional drift check
    if check_drift and not no_check_drift:
        result = drift_runner.run_drift_check(dataset)
        click.echo(f"drift: {result.status}")
        if result.status == "drifted" and result.diff:
            click.echo(json.dumps(result.diff, indent=2, default=str))
