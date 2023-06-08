"""Utilities for querying different databases for Target """
import lightkurve as lk
import astropy.units as u
from functools import lru_cache

@lru_cache
def get_timeseries(ra: u.Quantity, dec: u.Quantity) -> lk.LightCurve:
    """Function returns all the possible time series
     of an object as a Lightkurve object"""
    
    # query MAST for Kepler/TESS/K2

    # in theory we could grab WASP? ASAS-SN? all sorts

    #return lc
    raise NotImplementedError