__version__ = "0.1.0"
# Standard library
import os  # noqa

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

# Standard library
import logging  # noqa: E402

logging.basicConfig()
logger = logging.getLogger("exoscraper")
