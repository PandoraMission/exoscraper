# First-party/Local
import astropy.units as u

from exoscraper.planet import Planet, Planets


def test_planet():
    p = Planet("HAT-P-11", "b")
    p.PlanetParametersTable
    p.StarParametersTable
    assert isinstance(p.pl_rade, u.Quantity)
    assert hasattr(p.pl_rade, "err1")


def test_planets():
    p = Planets("HAT-P-11")
    p.PlanetsParametersTable
    p.StarParametersTable
    assert len(p) == 2
    assert isinstance(p.pl_rade, list)
