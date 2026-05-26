"""hvantk download <plugin> — top-level download commands.

Each plugin's downloader CLI is wired from its plugin.yaml ``cli:`` block
by the plugin loader. Adding a new plugin downloader does not require
editing this file — declare the entry under ``cli:`` in the plugin's
manifest (see hvantk/skills/_conventions/SKILL.md).

The loader strips the ``-download`` suffix from each ``command`` name to
derive the subcommand name within this group (e.g. ``clinvar-download``
becomes the ``clinvar`` subcommand under ``hvantk download``).
"""
import click

from hvantk.core.config import CONTEXT_SETTINGS
from hvantk.core.plugin import loader as plugin_loader


@click.group("download", context_settings=CONTEXT_SETTINGS)
def download_group():
    """Download external datasets."""


# Manifest-driven wiring: each plugin.yaml's cli: block declares its downloader.
# The loader strips the "-download" suffix from command names to derive the
# subcommand name within this group.
plugin_loader.get_registry().apply_plugin_downloaders(download_group)
