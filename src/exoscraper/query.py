"""Utilities for querying different databases for Target """
import lightkurve as lk
import astropy.units as u
from functools import lru_cache

@lru_cache
def get_timeseries(ra: u.Quantity, dec: u.Quantity) -> lk.LightCurve:
    """Function returns all the possible time series
     of an object as a Lightkurve object"""
    
    # query MAST for Kepler/TESS/K2

    # in theory we could grab WASP? ASAS-SN? ZTF? all sorts

    #return lc
    raise NotImplementedError

@lru_cache
def get_names(ra: u.Quantity, dec: u.Quantity) -> list:
    """Function to parse and retrieve all available names for a single target from Simbad"""

    # query simbad catalogs for ra and dec

    # return list of strings? There's gotta be a better format
    raise NotImplementedError

@lru_cache
def get_bibliography(names: list) -> dict: #?
    """Function to query NASA ADS for publications about this planet"""

    # parse names if names doesn't exist?
    # query NASA ADS based on names

    # return dictionary of references and links
    raise NotImplementedError

@lru_cache
def get_params(
        ra: u.Quantity,
        dec: u.Quantity,
        names: list,
        boundaries: dict,
) -> pd.DataFrame:
    """Function to query NASA Exoplanet Archive for planet parameters"""

    # query Exoplanet Archive for a set of parameters
    # if ra & dec are specified, fetch best match for those coords
    # same goes for names
    # if boundaries dict is specified, use those values to slice param space
    # perform some data validation to remove NaNs and unphysical values

    # return dictionary of parameters and values
    raise NotImplementedError

