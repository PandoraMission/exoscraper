"""Utilities for querying different databases for Target """
import warnings
from functools import lru_cache
from typing import List, Optional, Union

import astropy.units as u
import lightkurve as lk
import numpy as np
import pandas as pd
from astropy.io import votable
from astropy.constants import c as speedoflight
from astropy.utils.data import download_file
from astropy.coordinates import Distance, SkyCoord
from astropy.time import Time
from astroquery import log as asqlog
from astroquery.gaia import Gaia
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

asqlog.setLevel("ERROR")


@lru_cache
def get_SED(coord: Union[str, tuple], radius: Union[float, u.Quantity] = 2) -> dict:
    """Get the SED data for the target from Vizier

    Parameters
    ----------
    coord: string
        Astropy tuple of ra and dec or name of the object to query
    radius: float
        Radius to query in arcseconds
    """

    if isinstance(radius, u.Quantity):
        radius = radius.to(u.arcsecond).value
    if isinstance(coord, str):
        vizier_url = f"https://vizier.cds.unistra.fr/viz-bin/sed?-c={coord.replace(' ', '%20')}&-c.rs={radius}"
    elif isinstance(coord, tuple):
        vizier_url = f"https://vizier.cds.unistra.fr/viz-bin/sed?-c={coord[0]},{coord[1]}&-c.rs={radius}"
    else:
        raise ValueError("`coord` must be a `string` or `tuple` object.")

    df = votable.parse(download_file(vizier_url)).get_first_table().to_table()
    df = df[df["sed_flux"] / df["sed_eflux"] > 3]
    if len(df) == 0:
        # TODO: Need a logging loud warning here that there is no SED photometry
        return None
    wavelength = (speedoflight / (np.asarray(df["sed_freq"]) * u.GHz)).to(u.angstrom)
    sed_flux = np.asarray(df["sed_flux"]) * u.jansky
    sed_flux = sed_flux.to(
        u.erg / u.cm**2 / u.s / u.angstrom,
        equivalencies=u.spectral_density(wavelength),
    )
    sed_flux_err = np.asarray(df["sed_eflux"]) * u.jansky
    sed_flux_err = sed_flux_err.to(
        u.erg / u.cm**2 / u.s / u.angstrom,
        equivalencies=u.spectral_density(wavelength),
    )
    s = np.argsort(wavelength)
    SED = {
        "wavelength": wavelength[s],
        "sed_flux": sed_flux[s],
        "sed_flux_err": sed_flux_err[s],
        "filter": np.asarray(df["sed_filter"])[s],
    }
    return SED


@lru_cache
def get_timeseries(ra: u.Quantity, dec: u.Quantity) -> lk.LightCurve:
    """Function returns all the possible time series
    of an object as a Lightkurve object"""

    # query MAST for Kepler/TESS/K2

    # in theory we could grab WASP? ASAS-SN? ZTF? all sorts

    # return lc
    raise NotImplementedError


@lru_cache
def get_alternate_names(ra: u.Quantity, dec: u.Quantity) -> list:
    """Function to parse and retrieve all available names for a single target from Simbad"""

    # query simbad catalogs for ra and dec

    # return list of strings? There's gotta be a better format
    raise NotImplementedError


@lru_cache
def get_bibliography(names: list) -> dict:  # ?
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


@lru_cache
def get_sky_catalog(
    ra: float = 210.8023,
    dec: float = 54.349,
    radius: float = 0.155,
    gbpmagnitude_range: tuple = (-3, 20),
    limit: Optional = None,
    gaia_keys: List = [
        "source_id",
        "ra",
        "dec",
        "parallax",
        "pmra",
        "pmdec",
        "radial_velocity",
        "ruwe",
        "phot_bp_mean_mag",
        "teff_gspphot",
        "logg_gspphot",
    ],
) -> SkyCoord:
    """Gets a catalog of coordinates on the sky based on an input ra, dec and radius in units of degrees.

    RA and Dec are assumed to be in epoch J2000, this function will propagate the Gaia coordinates to J2000.
    """

    query_str = f"""
    SELECT {f'TOP {limit} ' if limit is not None else ''}* FROM (
        SELECT gaia.{', gaia.'.join(gaia_keys)}, dr2.teff_val AS dr2_teff_val,
        dr2.rv_template_logg AS dr2_logg, tmass.j_m, tmass.j_msigcom, tmass.ph_qual, DISTANCE(
        POINT({u.Quantity(ra, u.deg).value}, {u.Quantity(dec, u.deg).value}),
        POINT(gaia.ra, gaia.dec)) AS ang_sep,
        EPOCH_PROP_POS(gaia.ra, gaia.dec, gaia.parallax, gaia.pmra, gaia.pmdec,
        gaia.radial_velocity, gaia.ref_epoch, 2000) AS propagated_position_vector
        FROM gaiadr3.gaia_source AS gaia
        JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xmatch USING (source_id)
        JOIN gaiadr3.dr2_neighbourhood AS xmatch2 ON gaia.source_id = xmatch2.dr3_source_id
        JOIN gaiadr2.gaia_source AS dr2 ON xmatch2.dr2_source_id = dr2.source_id
        JOIN gaiadr3.tmass_psc_xsc_join AS xjoin USING (clean_tmass_psc_xsc_oid)
        JOIN gaiadr1.tmass_original_valid AS tmass ON
        xjoin.original_psc_source_id = tmass.designation
        WHERE 1 = CONTAINS(
        POINT({u.Quantity(ra, u.deg).value}, {u.Quantity(dec, u.deg).value}),
        CIRCLE(gaia.ra, gaia.dec, {(u.Quantity(radius, u.deg) + 50*u.arcsecond).value}))
        AND gaia.parallax IS NOT NULL
        AND gaia.phot_bp_mean_mag > {gbpmagnitude_range[0]}
        AND gaia.phot_bp_mean_mag < {gbpmagnitude_range[1]}) AS subquery
    WHERE 1 = CONTAINS(
    POINT({u.Quantity(ra, u.deg).value}, {u.Quantity(dec, u.deg).value}),
    CIRCLE(COORD1(subquery.propagated_position_vector), COORD2(subquery.propagated_position_vector), {u.Quantity(radius, u.deg).value}))
    ORDER BY ang_sep ASC
    """
    job = Gaia.launch_job_async(query_str, verbose=False)
    tbl = job.get_results()
    if len(tbl) == 0:
        raise ValueError("Could not find matches.")
    plx = tbl["parallax"].value.filled(fill_value=0)
    plx[plx < 0] = 0
    cat = {
        "jmag": tbl["j_m"].data.filled(np.nan),
        "bmag": tbl["phot_bp_mean_mag"].data.filled(np.nan),
        "ang_sep": tbl["ang_sep"].data.filled(np.nan) * u.deg,
    }
    cat["teff"] = (
        tbl["teff_gspphot"].data.filled(tbl["dr2_teff_val"].data.filled(np.nan)) * u.K
    )
    cat["logg"] = tbl["logg_gspphot"].data.filled(tbl["dr2_logg"].data.filled(np.nan))
    cat["RUWE"] = tbl["ruwe"].data.filled(99)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Double check the obstime is J2000
        cat["coords"] = SkyCoord(
            ra=tbl["ra"].value.data * u.deg,
            dec=tbl["dec"].value.data * u.deg,
            pm_ra_cosdec=tbl["pmra"].value.filled(fill_value=0) * u.mas / u.year,
            pm_dec=tbl["pmdec"].value.filled(fill_value=0) * u.mas / u.year,
            obstime=Time.strptime("2016", "%Y"),
            distance=Distance(parallax=plx * u.mas, allow_negative=True),
            radial_velocity=tbl["radial_velocity"].value.filled(fill_value=0)
            * u.km
            / u.s,
        ).apply_space_motion(Time.now())
    cat["source_id"] = np.asarray(
        [f"Gaia DR3 {i}" for i in tbl["source_id"].value.data]
    )
    return cat


@lru_cache
def get_planets(
    #    coord: SkyCoord,
    ra: float,
    dec: float,
    radius: u.Quantity = 20 * u.arcsecond,
    attrs: List = ["pl_orbper", "pl_tranmid", "pl_trandur", "pl_trandep"],
) -> dict:
    """
    Returns a dictionary of dictionaries with planet parameters.

    We assume RA and Dec are in J2000 epoch
    Largish default radius for high proper motion targets this breaks
    """
    # try:
    #     coord2000 = coord.apply_space_motion(Time(2000, format="jyear"))
    # except ValueError:
    #     coord2000 = coord
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        planets_tab = NasaExoplanetArchive.query_region(
            table="pscomppars", coordinates=SkyCoord(ra, dec, unit=u.deg), radius=radius
        )
        if len(planets_tab) != 0:
            planets = {
                letter: {
                    attr: planets_tab[planets_tab["pl_letter"] == letter][attr][
                        0
                    ].unmasked
                    for attr in attrs
                }
                for letter in planets_tab["pl_letter"]
            }
            # There's an error in the NASA exoplanet archive units that makes duration "days" instead of "hours"
            for planet in planets:
                if "pl_trandur" in planets[planet].keys():
                    planets[planet]["pl_trandur"] = (
                        planets[planet]["pl_trandur"].value * u.hour
                    )
        else:
            planets = {}
    return planets
