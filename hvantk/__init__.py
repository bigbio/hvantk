import logging

logger = logging.getLogger(__name__)
logger.info("Initializing hvantk package")


def main(*args, **kwargs):
    """Entrypoint proxy for `python -m hvantk`."""
    from .hvantk import main as _hvantk_main

    return _hvantk_main(*args, **kwargs)
