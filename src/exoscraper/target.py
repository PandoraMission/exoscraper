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
