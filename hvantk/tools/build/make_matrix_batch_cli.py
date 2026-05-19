"""
Batch builder for MatrixTables from a recipe file.

Recipe format (JSON or YAML):
{
  "matrices": [
    {
      "name": "ucsc-cellbrowser:default",
      "inputs": {"expression_matrix": "/path/expr.tsv.bgz", "metadata": "/path/meta.tsv"},
      "output": "/out/ucsc.mt",
      "params": {"gene_column": "gene", "overwrite": true}
    },
    {
      "name": "expression-atlas:dataset",
      "inputs": {"expression_matrix": "/path/matrix.tsv", "sdrf": "/path/atlas.sdrf.tsv"},
      "output": "/out/atlas.mt",
      "params": {"gene_column": "Gene ID"}
    }
  ]
}
"""

import json
import logging
from typing import Any, Dict

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)

try:
    import yaml  # optional
except Exception:  # pragma: no cover
    yaml = None


def _load_recipe(recipe_path: str) -> Dict[str, Any]:
    with open(recipe_path, "r") as fh:
        text = fh.read()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass
    if yaml is not None:
        return yaml.safe_load(text)
    raise ValueError(
        "Recipe is not valid JSON, and YAML support is not installed. Install PyYAML or provide JSON."
    )


@click.command("mkmatrix-batch", short_help="Create MatrixTables from a recipe.")
@click.option(
    "--recipe",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="JSON (or YAML if PyYAML installed) recipe describing matrices to build",
)
@click.pass_context
def mkmatrix_batch_cli(ctx, recipe):
    from hvantk.core.plugin.registry import run_matrix_builder

    try:
        spec = _load_recipe(recipe)
    except Exception as e:
        click.echo(str(e))
        ctx.abort()

    matrices = spec.get("matrices", [])
    if not isinstance(matrices, list) or not matrices:
        click.echo("Recipe is missing non-empty 'matrices' list")
        ctx.abort()

    for entry in matrices:
        name = entry.get("name")
        inputs = entry.get("inputs", {})
        output_mt = entry.get("output")
        params = entry.get("params", {})
        if not name or not isinstance(inputs, dict) or not output_mt:
            click.echo(
                f"Invalid entry in recipe (requires name,inputs,output): {entry}"
            )
            ctx.abort()
        logger.info(f"Building matrix via recipe: {name}")
        try:
            run_matrix_builder(name, inputs, output_mt, params)
        except Exception as e:
            click.echo(str(e))
            ctx.abort()
        click.echo(f"{name} MatrixTable created at {output_mt}")
