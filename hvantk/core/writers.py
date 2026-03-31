"""
Hail Table writer — universal persistence layer.

All algorithm outputs are written as Hail Tables regardless of the
compute backend used, ensuring downstream stages can always consume
them via any reader.
"""

import logging
from typing import List, Optional, Union

import pandas as pd

logger = logging.getLogger(__name__)


class HailTableWriter:
    """Write data as a Hail Table.

    Accepts either a pandas DataFrame or a Hail Table.  DataFrames are
    converted via ``hl.Table.from_pandas()`` before writing.
    """

    def write(
        self,
        data: Union[pd.DataFrame, "hl.Table"],
        path: str,
        key: Optional[List[str]] = None,
        overwrite: bool = False,
    ) -> None:
        """Persist data as a Hail Table.

        Parameters
        ----------
        data : pd.DataFrame or hl.Table
            The data to write.
        path : str
            Output ``.ht`` directory path.
        key : list[str], optional
            Fields to key the Hail Table by.
        overwrite : bool
            Overwrite an existing table at *path*.
        """
        import hail as hl

        if isinstance(data, pd.DataFrame):
            logger.info("Converting DataFrame (%d rows) to Hail Table", len(data))
            ht = hl.Table.from_pandas(data)
        else:
            ht = data

        if key:
            ht = ht.key_by(*key)

        logger.info("Writing Hail Table to %s", path)
        ht.write(path, overwrite=overwrite)
