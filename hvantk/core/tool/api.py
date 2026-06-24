"""Stable contracts for hvantk tools (analytical commands).

Tools are CLI commands that operate on data already inside hvantk's data
models. They are the analysis/transformation counterpart to plugins (which
bring data in). Each tool's tool.yaml is purely metadata; command registration
still happens via the existing top-level CLI wiring in hvantk/hvantk.py.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal


Domain = Literal[
    "plugins",      # plugin-system runtime
    "expression",
    "annotation",
    "genesets",
    "ptm",
    "ancestry",
    "qtl",
    "enrichex",
    "hgc",
    "infra",
    "rerank",
]

ToolType = Literal["command", "command_group"]


class ToolLoadError(Exception):
    """Raised when a tool.yaml fails to load, parse, or validate."""


@dataclass(frozen=True)
class Subcommand:
    name: str
    description: str
    inputs: tuple = ()                        # tuple of dicts {name, type, required, description}
    outputs: tuple = ()                       # tuple of output descriptors (strings or dicts)


@dataclass(frozen=True)
class ToolRequirements:
    hail: bool = False
    network: bool = False
    extras: tuple[str, ...] = ()              # additional Python package extras


@dataclass(frozen=True)
class ToolSpec:
    """One tool, materialised from a tool.yaml manifest."""

    name: str                                 # CLI verb (e.g., "drift", "plugins", "reprocess")
    domain: Domain
    type: ToolType
    description: str
    cli_module: str                            # e.g., "hvantk.tools.plugins.drift_cli"
    cli_callable: str                          # name of the Click group/command (e.g., "drift_cmd")
    purpose_short: str
    purpose_long: str = ""
    subcommands: tuple[Subcommand, ...] = ()
    requires: ToolRequirements = field(default_factory=ToolRequirements)
    manifest_path: str = ""                    # absolute path to the tool.yaml that produced this
