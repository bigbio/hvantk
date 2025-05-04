import logging

logger = logging.getLogger(__name__)
logger.info("Initializing hvantk package")

from .hvantk import main
