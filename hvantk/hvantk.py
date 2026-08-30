import importlib
import logging

import click

logger = logging.getLogger(__name__)

from hvantk.core.config import CONTEXT_SETTINGS

# Subcommand registry: name -> (module, attribute, short help).
#
# Importing the 18 subcommand modules eagerly pulled Hail, pandas, scipy and
# matplotlib into *every* invocation -- ``hvantk --help`` took ~11 s. They are
# resolved on first use instead (see ``LazyGroup``), so only the command the
# user actually ran pays for its dependencies.
#
# The short help is duplicated here because listing the commands must not import
# them. That duplication is the one maintenance cost of this design: if a
# subcommand's own help text changes, update the copy here too, or ``hvantk
# --help`` will show the stale wording.
_LAZY_COMMANDS: dict[str, tuple[str, str, str]] = {
    "ancestry-inference": (
        "hvantk.tools.ancestry.ancestry_cli",
        "ancestry_inference_cmd",
        "Infer genetic ancestry using PCA and Random Forest classification.",
    ),
    "annotate": (
        "hvantk.tools.annotation.annotate_cli",
        "annotate_group",
        "Build gene- and cohort-level annotation artifacts.",
    ),
    "catalog": (
        "hvantk.tools.infra.catalog_cli",
        "catalog",
        "Inspect the hvantk dataset catalog (aggregated from per-plugin catalogs).",
    ),
    "cohort": (
        "hvantk.tools.cohort.cohort_cli",
        "cohort_group",
        "Validate and attach external cohorts.",
    ),
    "download": (
        "hvantk.tools.plugins.download_cli",
        "download_group",
        "Download external datasets.",
    ),
    "drift": (
        "hvantk.tools.plugins.drift_cli",
        "drift_cmd",
        "Compare a plugin's live drift-probe fingerprint against the expected file.",
    ),
    "enrichex": (
        "hvantk.tools.enrichex",
        "enrichex_group",
        "Gene set enrichment analysis commands.",
    ),
    "expression": (
        "hvantk.tools.expression.summarize_expression_cli",
        "expression_group",
        "Expression AnnData analysis commands.",
    ),
    "genesets": (
        "hvantk.tools.genesets.genesets_cli",
        "genesets_group",
        "Extract or prepare gene set collections.",
    ),
    "hgc": (
        "hvantk.tools.hgc",
        "hgc_group",
        "HGC (Hail-based Genotype Combiner) commands for joint genotyping workflows",
    ),
    "plugins": (
        "hvantk.tools.plugins.plugins_cli",
        "plugins_group",
        "Inspect the hvantk plugin registry.",
    ),
    "psroc": (
        "hvantk.tools.ptm.psroc_cli",
        "psroc_cmd",
        "PSROC: Prediction Score ROC Analysis",
    ),
    "ptm": (
        "hvantk.tools.ptm.ptm_cli",
        "ptm_group",
        "Post-translational modification variant classification commands.",
    ),
    "qtlcascade": (
        "hvantk.tools.qtl.qtlcascade_cli",
        "qtlcascade_group",
        "Molecular QTL cascade analysis (eQTL \u2192 pQTL \u2192 disease).",
    ),
    "reprocess": (
        "hvantk.tools.plugins.reprocess_cli",
        "reprocess_cmd",
        "Run download -> parse -> build for a plugin dataset.",
    ),
    "rerank": (
        "hvantk.tools.rerank",
        "rerank_cmd",
        "Re-rank genes by multi-omic credibility from a declarative YAML config.",
    ),
    "tools": (
        "hvantk.tools.plugins.tools_cli",
        "tools_group",
        "Inspect the hvantk tool registry.",
    ),
    "utils": (
        "hvantk.tools.infra.utils_cli",
        "utils_group",
        "Operational utilities: format conversion, validation, diagnostics.",
    ),
}


class LazyGroup(click.Group):
    """A ``click.Group`` that imports a subcommand only when it is invoked.

    ``list_commands`` and ``--help`` are served from ``_LAZY_COMMANDS`` without
    importing anything; ``get_command`` does the real import and caches the
    resulting command object.
    """

    def list_commands(self, ctx):
        return sorted({*super().list_commands(ctx), *_LAZY_COMMANDS})

    def get_command(self, ctx, cmd_name):
        cmd = super().get_command(ctx, cmd_name)
        if cmd is not None:
            return cmd
        entry = _LAZY_COMMANDS.get(cmd_name)
        if entry is None:
            return None
        module, attr, _ = entry
        cmd = getattr(importlib.import_module(module), attr)
        self.add_command(cmd, cmd_name)
        return cmd

    def format_commands(self, ctx, formatter):
        """Render the command list from the registry, without importing."""
        rows = []
        for name in self.list_commands(ctx):
            entry = _LAZY_COMMANDS.get(name)
            if entry is not None:
                rows.append((name, entry[2]))
                continue
            cmd = super().get_command(ctx, name)
            if cmd is None or cmd.hidden:
                continue
            rows.append((name, cmd.get_short_help_str(limit=90)))
        if rows:
            with formatter.section("Commands"):
                formatter.write_dl(rows)


# Main CLI entry point for the package (hvantk)


def setup_logging(verbosity: int = 0, log_file: str | None = None):
    """Configure centralized logging for the hvantk CLI.

    Parameters
    ----------
    verbosity : int
        Verbosity level: 0 = WARNING (default), 1 = INFO, 2+ = DEBUG.
    log_file : str or None
        Optional path to a file where log output will be written.

    """
    level = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}.get(
        verbosity, logging.DEBUG
    )

    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=level,
        format="%(asctime)s %(name)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
        force=True,
    )


@click.group(
    "hvantk",
    cls=LazyGroup,
    help="A python package for gene and variant annotation with joint genotyping capabilities.",
    context_settings=CONTEXT_SETTINGS,
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase verbosity (-v: INFO, -vv: DEBUG)",
)
@click.option(
    "--log-file",
    type=click.Path(),
    default=None,
    help="Write logs to file",
)
def cli(verbose, log_file):
    """
    Entry point for the hvantk command-line interface.

    Serves as the root CLI group for gene and variant annotation commands
    with integrated joint genotyping workflows.
    """
    setup_logging(verbose, log_file)
    logger.info("Starting hvantk CLI")


def main():
    """
    Runs the main entry point for the hvantk CLI application.

    Invokes the command-line interface to process user commands and logs the start and completion of the main function.
    """
    logger.info("Running main function")
    cli()
    logger.info("Main function completed")


if __name__ == "__main__":
    main()
