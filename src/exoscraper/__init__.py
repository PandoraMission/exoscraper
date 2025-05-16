__version__ = "0.2.1"
# Standard library
import os  # noqa

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

# Standard library
import logging  # noqa: E402

# This library lets us have log messages with syntax highlighting
from rich.logging import RichHandler  # noqa: E402

# logging.basicConfig()
# logger = logging.getLogger("exoscraper")
log = logging.getLogger("exoscraper")
log.addHandler(RichHandler(markup=True))
log.setLevel("INFO")

from .system import System  # noqa
