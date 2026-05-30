"""Plugin-local constants for the alphagenome skill.

Moved from hvantk/core/constants.py per issue #120 (clean-core principle).
"""

ALPHAGENOME_DEFAULT_INTERVAL_SIZE = 1_048_576  # 1Mbp
ALPHAGENOME_DEFAULT_DENSITY_WINDOW = 50_000    # 50kb
ALPHAGENOME_DEFAULT_RETRY_BACKOFF = 2.0
ALPHAGENOME_DEFAULT_MAX_RETRIES = 3
ALPHAGENOME_DEFAULT_REQUEST_TIMEOUT = 120
