"""Classes for working with Targets"""
from dataclasses import dataclass
from typing import Union, List
from astropy.coordinate import SkyCoord
from query import get_timeseries

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
        return

    @staticmethod
    def from_gaia(coord: Union[str, SkyCoord]):
        raise NotImplementedError
    
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
    def from_names(coords: List[str | SkyCoord]):
        raise NotImplementedError
    
    def to_csv(self, output):
        """Produces csv file with all the targets in the TargetSet """
        raise NotImplementedError
