# First-party/Local
from exoscraper import Target, __version__


def test_version():
    assert __version__ == "0.1.0"


def test_init():
    """Check that Target can be initialized"""
    t = Target.from_gaia("HAT-P-19")
    assert isinstance(t.planets, dict)


#    assert ('b' is in t.planets.keys)
