import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


@click.group("utils", context_settings=CONTEXT_SETTINGS)
def utils_group():
    """Operational utilities: format conversion, validation, diagnostics."""


from hvantk.tools.infra.check_install_cli import check_install_cmd
from hvantk.tools.build.build_1k_genome_cli import build_1k_genome_cmd
from hvantk.tools.infra.validate_bgzf_cli import validate_bgzf_cmd

utils_group.add_command(check_install_cmd)
utils_group.add_command(build_1k_genome_cmd)
utils_group.add_command(validate_bgzf_cmd)

try:
    from hvantk.tools.infra.convert_bgz_cli import convert_bgz_cmd

    utils_group.add_command(convert_bgz_cmd)
except ImportError as e:
    logger.debug(
        "Could not load convert-bgz command (%s): %s",
        "hvantk.tools.infra.convert_bgz_cli",
        e,
    )
