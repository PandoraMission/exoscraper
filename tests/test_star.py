# First-party/Local
import astropy.units as u

from exoscraper.star import Star
from exoscraper.query import get_planets


def test_star():
    """Check to see if Star can be initialized and values read out"""
    tab = get_planets(name="HAT-P-19")
    s = Star(tab)
    assert isinstance(s.st_rad, u.Quantity)
    assert hasattr(s.st_rad, "err1")
    assert len(s.planets) >= 1
