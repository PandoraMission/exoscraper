# First-party/Local
from exoscraper import system, __version__

import numpy as np


def test_version():
    assert __version__ == "0.2.1"


def test_init():
    """Check that Target can be initialized"""
    t = system.System.from_gaia("HAT-P-19")
    assert isinstance(t.planets, list)
    t = system.System.from_name("HAT-P-19")
    assert isinstance(t.planets, list)
    assert "HAT-P-19b" in t.planets[0].name

    # Can get SED property
    sed = t.SED
    assert isinstance(sed, dict)
    assert len(sed["wavelength"]) > 0


def test_timeseries_model():
    """Check that the BATMAN implementation to model a time series works"""
    t = system.System.from_gaia("TOI-700")
    time = np.linspace(2458000, 2458100, 10000)
    model0, vars0 = t.sample_timeseries(
        time=time, median_flag=True, vals_out=True, iterations=3
    )
    assert isinstance(vars0, dict)
    assert isinstance(model0, np.ndarray)
    assert model0.shape == (3, 10000)

    # Check that overwriting median values with user-provided ones works
    model, vars = t.sample_timeseries(
        time=time, median_flag=True, vals_out=True, iterations=3, period=[1, 2, 4, 8]
    )
    assert sum(sum(model0 == model)) != 30000
    assert (vars0["period"] != vars["period"]).all()
