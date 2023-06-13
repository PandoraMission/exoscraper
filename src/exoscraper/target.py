"""Classes for working with Targets"""
from dataclasses import dataclass
from typing import Union

import astropy.units as u
import lightkurve as lk
from astropy.coordinates import SkyCoord

from .query import get_planets, get_sky_catalog


@dataclass
class Target(object):
    name: str
    ra: u.Quantity
    dec: u.Quantity
    logg: u.Quantity
    teff: u.Quantity
    bmag: u.Quantity
    jmag: u.Quantity
    coord: SkyCoord = None
    """DOCSTRING"""

    def __post_init__(self):
        """Ensures quantity conventions"""
        self.ra, self.dec = u.Quantity(self.ra, u.deg), u.Quantity(self.dec, u.deg)
        self.teff = u.Quantity(self.teff, u.K)
        self.logg = u.Quantity(self.logg)
        self.planets = get_planets(self.ra.value, self.dec.value)
        return

    def __repr__(self):
        return f"{self.name} [{self.ra}, {self.dec}]"

    def _repr_html_(self):
        return f"{self.name} ({self.ra._repr_latex_()},  {self.dec._repr_latex_()})"

    @staticmethod
    def from_gaia(coord: Union[str, SkyCoord]):
        name = None
        if isinstance(coord, str):
            name = coord
            coord = SkyCoord.from_name(coord)
        elif not isinstance(coord, SkyCoord):
            raise ValueError("`coord` must be a `SkyCoord` or a name string.")
        cat = get_sky_catalog(coord.ra, coord.dec, radius=5 * u.arcsecond, limit=1)
        if name is None:
            name = cat["source_id"][0]
        return Target(
            name=name,
            ra=cat["coords"][0].ra,
            dec=cat["coords"][0].dec,
            logg=cat["logg"][0],
            teff=cat["teff"][0],
            bmag=cat["bmag"][0],
            jmag=cat["jmag"][0],
            coord=cat["coords"][0],
        )

    @staticmethod
    def from_TIC(coord: Union[str, SkyCoord]):
        raise NotImplementedError

    @staticmethod
    def from_name(coord: str):
        raise NotImplementedError

    @property
    def lightcurve(self):
        # go get the TESS data something like
        # return get_timeseries(self.ra, self.dec)
        raise NotImplementedError

    @property
    def bibliography(self):
        # go get references from NASA ADS
        # return get_bibliography(self.name)
        raise NotImplementedError

    @property
    def transit_times(self, time: u.Quantity):
        # calculate future transit times
        # return array of future transit times out to a specified time
        raise NotImplementedError

    @property
    def SED(self):
        # calculate SED of target star
        # return idk whatever Christina has written
        raise NotImplementedError

    @property
    def transit_model(self) -> lk.LightCurve:
        # generate transit model using lightkurve
        # return LightCurve object
        raise NotImplementedError

    @property
    def noise_model(self, mission: str) -> float:
        # generate noise model for target for specified mission
        # this could be wrapped into transit model
        # this also may not be worth including since I selfishly need it for
        #  pandora-target and so added it here lol
        # returns LightCurve object
        raise NotImplementedError

    def motto(self):
        # mottos = ["If it's online, we can scrape it!",
        #           "Leave the scraping to us!",
        #           "Scraping the internet clean!",
        #           "Scrape it or leave it!",
        #           "Consider it scraped!"]

        # print(random.choice(mottos))
        raise NotImplementedError


class TargetSet(object):
    """A class to hold many Target classes"""

    def __init__(self, targets):
        self.targets = targets

    def __iter__(self):
        raise NotImplementedError

    def __len__(self):
        raise NotImplementedError

    def __repr__(self):
        raise NotImplementedError

    @staticmethod
    def from_names(coords: Union[str, SkyCoord]):
        raise NotImplementedError

    def to_csv(self, output):
        """Produces csv file with all the targets in the TargetSet"""
        raise NotImplementedError
