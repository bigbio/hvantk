"""AlphaGenome variant prediction streamer.

Streams variant positions through the AlphaGenome API and produces
per-modality Hail Tables with full multimodal predictions.
"""

import json
import logging
import os
import random
import time
from dataclasses import dataclass
from typing import Any, Dict, Iterator, List, Optional, Tuple

import yaml

from hvantk.data.data_streamer import HailDataStreamer

logger = logging.getLogger(__name__)

_REQUIRED_SECTIONS = ("api", "ontology")
_DEFAULT_INTERVALS = {
    "default_size": 1_048_576,
    "adaptive": True,
    "adaptive_max_size": 1_048_576,
    "density_window": 50_000,
}


def load_config(config_path: str) -> Dict[str, Any]:
    """Load and validate an AlphaGenome YAML config file.

    Auth resolution order: config api.key > ALPHAGENOME_API_KEY env var > error.

    Parameters
    ----------
    config_path : str
        Path to YAML config file.

    Returns
    -------
    dict
        Validated config dict with resolved API key and interval defaults.

    Raises
    ------
    FileNotFoundError
        If config_path does not exist.
    ValueError
        If required sections are missing or API key cannot be resolved.
    """
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path) as f:
        config = yaml.safe_load(f)

    for section in _REQUIRED_SECTIONS:
        if section not in config or config[section] is None:
            raise ValueError(
                f"Config missing required section: '{section}'. "
                f"Required sections: {_REQUIRED_SECTIONS}"
            )

    # Resolve API key: config > env var > error
    api_key = config["api"].get("key")
    if not api_key:
        api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        raise ValueError(
            "API key not found. Set 'api.key' in config or "
            "ALPHAGENOME_API_KEY environment variable."
        )
    config["api"]["key"] = api_key

    # Apply interval defaults
    if "intervals" not in config or config["intervals"] is None:
        config["intervals"] = dict(_DEFAULT_INTERVALS)
    else:
        for k, v in _DEFAULT_INTERVALS.items():
            config["intervals"].setdefault(k, v)

    return config
