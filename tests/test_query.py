import astropy.units as u
from astropy.coordinates import SkyCoord

from exoscraper.query import get_SED


def test_sed():
    sed = get_SED("HAT-P-19", 1)
    assert isinstance(sed, dict)
    assert len(sed["wavelength"]) > 0
    assert isinstance(sed["wavelength"], u.Quantity)

    sed = get_SED("HAT-P-19", 1 * u.arcsecond)
    assert isinstance(sed, dict)
    assert len(sed["wavelength"]) > 0
    assert isinstance(sed["wavelength"], u.Quantity)

    coord = SkyCoord.from_name("HAT-P-19")
    sed = get_SED((coord.ra.value, coord.dec.value), 1)
    assert isinstance(sed, dict)
    assert len(sed["wavelength"]) > 0
    assert isinstance(sed["wavelength"], u.Quantity)


def test_bad_sed(caplog):
    # If we pass in an invalid name we should get no results
    sed = get_SED("christina", 1)
    assert sed is None
    assert "No SED photometry found" in caplog.text
