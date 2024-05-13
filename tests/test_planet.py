# First-party/Local
import astropy.units as u

from exoscraper.planet import Planet
from exoscraper.query import get_planets


def test_planet():
    """Check that a Planet instance can be initialized"""
    tab = get_planets(name="HAT-P-19")
    p = Planet(tab[0])
    p.PlanetParametersTable
    # p.StarParametersTable
    assert isinstance(p.pl_rade, u.Quantity)
    assert hasattr(p.pl_rade, "err1")


def test_distributions():
    """Check that distributions are working in the Planet class"""
    tab = get_planets(name="HAT-P-19")
    p = Planet(tab[0])
    assert len(p.pl_rade.distribution.sample(size=3)) == 3
    assert p.pl_rade.distribution.disttype == "Normal"


# def test_planets():
#     p = Planets("HAT-P-11")
#     p.PlanetsParametersTable
#     p.StarParametersTable
#     assert len(p) == 2
#     assert isinstance(p.pl_rade, list)
