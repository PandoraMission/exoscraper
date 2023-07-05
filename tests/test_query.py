from exoscraper.query import get_SED
import astropy.units as u
from astropy.coordinates import SkyCoord


def test_sed():
    sed = get_SED("HAT-P-19", 1)
    assert isinstance(sed, dict)
    assert len(sed["wavelength"]) > 0
    assert isinstance(sed["wavelength"], u.Quantity)

    coord = SkyCoord.from_name("HAT-P-19")
    sed = get_SED((coord.ra.value, coord.dec.value), 1)
    assert isinstance(sed, dict)
    assert len(sed["wavelength"]) > 0
    assert isinstance(sed["wavelength"], u.Quantity)
