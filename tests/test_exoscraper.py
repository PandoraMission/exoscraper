# First-party/Local
from exoscraper import Target, __version__


def test_version():
    assert __version__ == "0.1.0"


def test_init():
    """Check that Target can be initialized"""
    t = Target.from_gaia("HAT-P-19")
    assert isinstance(t.planets, dict)
    t = Target.from_name("HAT-P-19")
    assert isinstance(t.planets, dict)
    assert "b" in t.planets.keys()

    # Can get SED property
    sed = t.SED
    assert isinstance(sed, dict)
    assert len(sed["wavelength"]) > 0
