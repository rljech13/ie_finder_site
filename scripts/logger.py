"""Logging for finder pipeline scripts.

Messages go to stderr only. Snakemake redirects each step's output to that
step's log, and the report keeps any log that is not empty.
"""

from __future__ import annotations

import logging

FORMAT = "[%(asctime)s - %(levelname)s]: %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def get_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """Return a named logger that writes to stderr at ``level`` and above."""
    logging.basicConfig(format=FORMAT, datefmt=DATE_FORMAT)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    return logger
