"""Catalog CLI for inspecting and invoking artifact builders.

Supports listing entries, showing a single entry's metadata, and invoking its
registered builder callable (if available). The catalog YAML may be either:
- A list of entries (each with an 'id')  [current hvantk/resources/catalog.yaml]
- A mapping of id -> entry (legacy / alternate style)

Each entry may specify a builder at entry['provenance']['builder'] (current
format) or at entry['builder'] (alternate). Builder arguments are taken from
entry['args'] or entry['params'] (whichever present) and can be overridden via
--override key=value pairs.
"""

import importlib
import importlib.resources as ir
from typing import Any, Dict, List, Tuple
import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)

try:  # pragma: no cover - exercised in runtime environment
    import yaml  # type: ignore
except Exception:  # pragma: no cover
    yaml = None  # type: ignore


def _resolve_callable(dotted: str):
    mod, attr = dotted.rsplit(".", 1)
    m = importlib.import_module(mod)
    fn = getattr(m, attr, None)
    if not callable(fn):  # pragma: no cover - defensive
        raise RuntimeError(f"Builder {dotted} is not callable")
    return fn


def _normalize_catalog(raw: Any) -> Dict[str, Dict[str, Any]]:
    """Return a mapping id -> entry regardless of input shape."""
    if raw is None:
        return {}
    if isinstance(raw, dict) and all(isinstance(v, dict) for v in raw.values()):
        return raw  # already keyed
    if isinstance(raw, list):
        out: Dict[str, Dict[str, Any]] = {}
        for entry in raw:
            if not isinstance(entry, dict):  # pragma: no cover
                continue
            entry_id = entry.get("id")
            if not entry_id:
                continue
            out[entry_id] = entry
        return out
    # Fallback: unsupported shape
    logger.warning("Unrecognized catalog YAML structure; expected list or mapping")
    return {}


def _load_catalog(path: str | None = None) -> Dict[str, Dict[str, Any]]:
    if yaml is None:
        raise RuntimeError(
            "PyYAML is required for catalog CLI. Install with `pip install pyyaml`."
        )
    if path:
        with open(path, "r") as fh:
            raw = yaml.safe_load(fh)
    else:
        with ir.files("hvantk.resources").joinpath("catalog.yaml").open("r") as fh:
            raw = yaml.safe_load(fh)
    return _normalize_catalog(raw)


@click.group(
    "catalog",
    context_settings=CONTEXT_SETTINGS,
    help="Interact with hvantk artifact catalog (list, show, build).",
)
def catalog():
    """Root command group for catalog operations."""
    pass


@catalog.command("list")
@click.option(
    "--catalog",
    "catalog_path",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    default=None,
    help="Path to catalog YAML (defaults to packaged hvantk.resources.catalog.yaml)",
)
def list_entries(catalog_path: str | None):
    """List available catalog entry IDs."""
    cat = _load_catalog(catalog_path)
    if not cat:
        click.echo("(empty catalog)")
        return
    for k in sorted(cat.keys()):
        click.echo(k)


@catalog.command("show")
@click.argument("key")
@click.option(
    "--catalog",
    "catalog_path",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    default=None,
)
def show_entry(key: str, catalog_path: str | None):
    """Show full YAML for an entry."""
    cat = _load_catalog(catalog_path)
    entry = cat.get(key)
    if not entry:
        raise click.ClickException(f"Entry not found: {key}")
    click.echo(yaml.safe_dump(entry, sort_keys=False))


def _extract_builder(entry: Dict[str, Any]) -> Tuple[str, Dict[str, Any]]:
    # Prefer provenance.builder then top-level builder
    builder_path = entry.get("provenance", {}).get("builder") or entry.get("builder")
    if not builder_path:
        raise click.ClickException(
            "Entry has no builder path defined (provenance.builder or builder)"
        )
    # Args precedence: args > params
    args = dict(entry.get("args") or entry.get("params") or {})
    return builder_path, args


@catalog.command("build")
@click.argument("key")
@click.option(
    "--catalog",
    "catalog_path",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    default=None,
)
@click.option(
    "--override",
    multiple=True,
    help="Override builder args as key=value. Repeatable.",
)
def build_entry(key: str, catalog_path: str | None, override: List[str]):
    """Invoke the builder for a catalog entry (experimental)."""
    cat = _load_catalog(catalog_path)
    entry = cat.get(key)
    if not entry:
        raise click.ClickException(f"Entry not found: {key}")
    builder_path, args = _extract_builder(entry)
    for kv in override:
        if "=" not in kv:
            raise click.ClickException(f"Invalid override (expected key=value): {kv}")
        k, v = kv.split("=", 1)
        args[k] = v
    click.echo(f"Invoking builder {builder_path} with args: {args}")
    builder = _resolve_callable(builder_path)
    result = builder(**args)
    click.echo(f"Built {key}: {result}")
