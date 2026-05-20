"""HGC CLI Commands Module - Main Entry Point.

Exposes a Click group whose subcommands import the hail-dependent
implementation lazily. This lets ``hvantk --help`` and unrelated
subcommands work in environments where hail is not installed; the heavy
imports only fire when the user actually runs ``hvantk hgc ...``.
"""

import logging

import click

logger = logging.getLogger(__name__)


class _LazyHgcGroup(click.Group):
    """Click group that registers its subcommands on first access."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._loaded = False

    def _ensure_loaded(self) -> None:
        if self._loaded:
            return
        from .combine_cli import register_combine_commands
        from .convert_cli import register_convert_commands
        from .pipeline_cli import register_pipeline_command
        from .qc_cli import register_qc_commands

        register_combine_commands(self)
        register_convert_commands(self)
        register_qc_commands(self)
        register_pipeline_command(self)
        self._loaded = True

    def list_commands(self, ctx):  # type: ignore[override]
        self._ensure_loaded()
        return super().list_commands(ctx)

    def get_command(self, ctx, name):  # type: ignore[override]
        self._ensure_loaded()
        return super().get_command(ctx, name)


@click.group(
    name="hgc",
    cls=_LazyHgcGroup,
    help="HGC (Hail-based Genotype Combiner) commands for joint genotyping workflows",
)
@click.pass_context
def hgc_group(ctx):
    """HGC command group for joint genotyping operations.

    Logging is configured centrally via ``hvantk -v`` / ``hvantk --log-file``.
    """
    ctx.ensure_object(dict)
    logger.info("Starting HGC command")
