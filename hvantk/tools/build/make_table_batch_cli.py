"""
Generate annotation tables from a recipe (no legacy raw_data_path mode).
"""

import logging
from typing import Any, Dict

import click
import json

logger = logging.getLogger(__name__)

from hvantk.core.config import CONTEXT_SETTINGS

try:
    import yaml  # optional, only if installed
except Exception:  # pragma: no cover
    yaml = None


def _load_recipe(recipe_path: str) -> Dict[str, Any]:
    """Load a recipe file (JSON or YAML if PyYAML is installed)."""
    with open(recipe_path, "r") as fh:
        text = fh.read()
    # Try JSON first
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass
    # Try YAML if available
    if yaml is not None:
        return yaml.safe_load(text)
    raise ValueError(
        "Recipe file is not valid JSON, and YAML support is not installed. Install PyYAML or provide JSON."
    )


@click.command("mktable-batch", short_help="Create annotation tables from a recipe.")
@click.option(
    "--recipe",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="JSON (or YAML if PyYAML installed) recipe describing tables to build",
)
@click.pass_context
def mktable_batch_cli(ctx, recipe):
    """Create one or more annotation tables described in a recipe file."""
    from hvantk.core.plugin.registry import run_table_builder

    try:
        spec = _load_recipe(recipe)
    except Exception as e:
        click.echo(str(e))
        ctx.abort()

    tables = spec.get("tables", [])
    if not isinstance(tables, list) or not tables:
        click.echo("Recipe is missing non-empty 'tables' list")
        ctx.abort()

    for entry in tables:
        name = entry.get("name")
        input_path = entry.get("input")
        output_path = entry.get("output")
        params = entry.get("params", {})
        if not name or not input_path or not output_path:
            click.echo(f"Invalid entry in recipe (requires name,input,output): {entry}")
            ctx.abort()
        logger.info(f"Building table via recipe: {name}")
        run_table_builder(name, input_path, output_path, params)
        click.echo(f"{name} table created at {output_path}")


if __name__ == "__main__":
    mktable_batch_cli()  # pragma: no cover
