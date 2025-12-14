"""Utilities for querying different databases for Target"""

import warnings
from functools import lru_cache
from typing import List, Union
import json

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.constants import c as speedoflight
from astropy.coordinates import Distance, SkyCoord
from astropy.io import votable
from astropy.table import QTable
from astropy.time import Time
from astropy.utils.data import download_file
from astroquery import log as asqlog
from astroquery.gaia import Gaia
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from gaiaoffline import Gaia as Gaia_off

# from bs4 import BeautifulSoup
import requests

from . import log

# from .utils import calc_separation

asqlog.setLevel("ERROR")


@lru_cache
def get_SED(coord: Union[str, tuple], radius: Union[float, u.Quantity] = 2.0) -> dict:
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
    elif isinstance(radius, (int, float)):
        radius = float(radius)
    if isinstance(coord, str):
        vizier_url = f"https://vizier.cds.unistra.fr/viz-bin/sed?-c={coord.replace(' ', '%20')}&-c.rs={radius}"
    elif isinstance(coord, tuple):
        vizier_url = f"https://vizier.cds.unistra.fr/viz-bin/sed?-c={coord[0]},{coord[1]}&-c.rs={radius}"
    else:
        raise ValueError("`coord` must be a `string` or `tuple` object.")
    try:
        df = (
            votable.parse(download_file(vizier_url, show_progress=False))
            .get_first_table()
            .to_table()
        )
    except IndexError:
        log.warning(f"No SED photometry found for `{coord}` at Vizier.")
        return None

    df = df[df["sed_flux"] / df["sed_eflux"] > 3]
    if len(df) == 0:
        log.warning(f"No SED photometry found for {coord} at Vizier.")
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
def get_timeseries(ra: u.Quantity, dec: u.Quantity) -> np.ndarray:
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
    ra: float,
    dec: float,
    radius: float = 0.155,
    gbpmagnitude_range: tuple = (-3, 20),
    limit=None,
    gaia_keys: list = [],
    time: Time = Time.now(),
) -> dict:
    """
    Gets a catalog of coordinates on the sky based on an input RA, Dec, and radius as well as
    a magnitude range for Gaia. The user can also specify additional keywords to be grabbed
    from Gaia catalog.

    Parameters
    ----------
    ra : float
        Right Ascension of the center of the query radius in degrees.
    dec : float
        Declination of the center of the query radius in degrees.
    radius : float
        Radius centered on ra and dec that will be queried in degrees.
    gbpmagnitude_range : tuple
        Magnitude limits for the query. Targets outside of this range will not be included in
        the final output dictionary.
    limit : int
        Maximum number of targets from query that will be included in output dictionary. If a
        limit is specified, targets will be included based on proximity to specified ra and dec.
    gaia_keys : list
        List of additional Gaia archive columns to include in the final output dictionary.
    time : astropy.Time object
        Time at which to evaluate the positions of the targets in the output dictionary.

    Returns
    -------
    cat : dict
        Dictionary of values from the Gaia archive for each keyword.
    """

    base_keys = [
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
        "phot_g_mean_flux",
        "phot_g_mean_mag",
    ]

    all_keys = base_keys + gaia_keys

    query_str = f"""
    SELECT {f'TOP {limit} ' if limit is not None else ''}* FROM (
        SELECT gaia.{', gaia.'.join(all_keys)}, dr2.teff_val AS dr2_teff_val,
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
    job = Gaia.launch_job_async(query_str, verbose=True)
    tbl = job.get_results()
    if len(tbl) == 0:
        raise ValueError("Could not find matches.")
    plx = tbl["parallax"].value.filled(fill_value=0)
    plx[plx < 0] = 0
    cat = {
        "jmag": tbl["j_m"].data.filled(np.nan),
        "bmag": tbl["phot_bp_mean_mag"].data.filled(np.nan),
        "gmag": tbl["phot_g_mean_mag"].data.filled(np.nan),
        "gflux": tbl["phot_g_mean_flux"].data.filled(np.nan),
        "ang_sep": tbl["ang_sep"].data.filled(np.nan) * u.deg,
    }
    cat["teff"] = (
        tbl["teff_gspphot"].data.filled(tbl["dr2_teff_val"].data.filled(np.nan)) * u.K
    )
    cat["logg"] = tbl["logg_gspphot"].data.filled(tbl["dr2_logg"].data.filled(np.nan))
    cat["RUWE"] = tbl["ruwe"].data.filled(99)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cat["coords"] = SkyCoord(
            ra=tbl["ra"].value.data * u.deg,
            dec=tbl["dec"].value.data * u.deg,
            pm_ra_cosdec=tbl["pmra"].value.filled(fill_value=0) * u.mas / u.year,
            pm_dec=tbl["pmdec"].value.filled(fill_value=0) * u.mas / u.year,
            obstime=Time("2457389.0", format="jd", scale="tcb"),  # "J2016.0"
            distance=Distance(parallax=plx * u.mas, allow_negative=True),
            radial_velocity=tbl["radial_velocity"].value.filled(fill_value=0)
            * u.km
            / u.s,
        ).apply_space_motion(time)
    cat["source_id"] = np.asarray(
        [f"Gaia DR3 {i}" for i in tbl["source_id"].value.data]
    )
    for key in gaia_keys:
        cat[key] = tbl[key].data.filled(np.nan)
    return cat


def get_offline_star_catalog(
    target_ra: float,
    target_dec: float,
    radius: float = 0.155,
    input_epoch: float = 2000.0,
    limit=None,
    time: Time = Time.now(),
):
    """
    Gets a catalog of coordinates on the sky based on an input RA, Dec, and radius for Gaia
    DR3. This is meant to be a fully offline version of `get_sky_catalog` that accesses
    the Gaia DR3 catalog via the `gaiaoffline` package. Since it's offline, there are
    fewer available keys that can be accessed compared to the online query.

    Parameters
    ----------
    target_ra : float
        Right Ascension of the center of the query radius in degrees.
    target_dec : float
        Declination of the center of the query radius in degrees.
    radius : float
        Radius centered on ra and dec that will be queried in degrees.
    limit : int
        Maximum number of targets from query that will be included in output dictionary. If a
        limit is specified, targets will be included based on proximity to specified ra and dec.
    time : astropy.Time object
        Time at which to evaluate the positions of the targets in the output dictionary.

    Returns
    -------
    cat : dict
        Dictionary of values from the Gaia archive for each keyword.
    """
    # Do some input coord processing
    target_ra = u.Quantity(target_ra, u.deg).value
    target_dec = u.Quantity(target_dec, u.deg).value

    # Query the gaiaoffline DR3 catalog for targets around specified RA/Dec
    with Gaia_off(limit=1000, photometry_output="mag", tmass_crossmatch=True) as gaia:
        cat = gaia.conesearch(
            ra=target_ra,
            dec=target_dec,
            radius=radius,
        )

    cat_prop = cat.copy()

    # Propagate coordinates back to the epoch of the input coords
    # Handle missing radial velocities
    rv = cat_prop["radial_velocity"].fillna(0).values * u.km / u.s

    # Gaia DR3 reference epoch
    if "ref_epoch" in cat_prop.columns:
        ref_epoch = Time(cat_prop["ref_epoch"].values, format="jyear")
    else:
        ref_epoch = Time(2016.0, format="jyear")

    input_time = Time(input_epoch, format="jyear")

    # Turn DR3 entries into SkyCoord objects
    coords = SkyCoord(
        ra=cat_prop["ra"].values * u.deg,
        dec=cat_prop["dec"].values * u.deg,
        distance=Distance(
            parallax=cat_prop["parallax"].fillna(np.nan).values * u.mas,
            allow_negative=True,
        ),
        # distance=(1000.0 / cat_prop['parallax'].values) * u.pc,  # Convert parallax to distance
        pm_ra_cosdec=cat_prop["pmra"].fillna(0).values * u.mas / u.yr,
        pm_dec=cat_prop["pmdec"].fillna(0).values * u.mas / u.yr,
        radial_velocity=rv,
        obstime=ref_epoch,
        frame="icrs",
    )

    # Propagate to input epoch
    coords_input_propagated = coords.apply_space_motion(new_obstime=input_time)

    # Extract propagated coordinates
    cat_prop["input_propagated_ra"] = coords_input_propagated.ra.deg
    cat_prop["input_propagated_dec"] = coords_input_propagated.dec.deg

    # Create propagated coordinates
    prop_coords = SkyCoord(
        ra=cat_prop["input_propagated_ra"].values * u.deg,
        dec=cat_prop["input_propagated_dec"].values * u.deg,
        frame="icrs",
    )

    # Create target coordinate
    target_coord = SkyCoord(ra=target_ra * u.deg, dec=target_dec * u.deg, frame="icrs")

    # Determine the coordinate set nearest the input at the input epoch
    # Calculate angular separations at propagated epoch
    separations = target_coord.separation(prop_coords)
    cat_prop["input_epoch_separation_arcsec"] = separations.arcsec

    # Filter by search radius at propagated epoch
    mask = cat_prop["input_epoch_separation_arcsec"] <= 5.0
    matches = cat_prop[mask]

    if len(matches) == 0:
        print(f"No stars found within {5.0} arcsec at epoch {input_epoch}")
        return None

    # Find the nearest match
    nearest_idx = matches["input_epoch_separation_arcsec"].idxmin()

    # Propagate the coordinates to the specified input time
    new_coords = coords.apply_space_motion(time)
    cat_prop["new_coords"] = new_coords

    # Create updated target coordinate based on matching
    target_coord = new_coords[nearest_idx]

    separations = target_coord.separation(new_coords)
    cat_prop["output_epoch_separation_arcsec"] = separations.arcsec
    cat_prop = cat_prop.sort_values(by="output_epoch_separation_arcsec").reset_index(
        drop=True
    )

    # Build the output catalog
    out_cat = {
        "jmag": cat_prop["j_m"].fillna(np.nan),
        "bmag": cat_prop["phot_bp_mean_mag"].fillna(np.nan),
        "gmag": cat_prop["phot_g_mean_mag"].fillna(np.nan),
        # "gflux": cat_prop["phot_g_mean_flux"].fillna(np.nan),
    }
    out_cat["ang_sep"] = cat_prop["output_epoch_separation_arcsec"]
    out_cat["teff"] = cat_prop["teff_gspphot"].fillna(np.nan).tolist() * u.K
    out_cat["logg"] = cat_prop["logg_gspphot"].fillna(np.nan).tolist()
    # out_cat["RUWE"] = cat_prop["ruwe"].fillna(99)

    ras = [coord.ra.deg for coord in cat_prop["new_coords"]]
    decs = [coord.dec.deg for coord in cat_prop["new_coords"]]
    distances = [coord.distance.pc for coord in cat_prop["new_coords"]]
    pm_ra_cosdec = [coord.pm_ra_cosdec.value for coord in cat_prop["new_coords"]]
    pm_dec = [coord.pm_dec.value for coord in cat_prop["new_coords"]]
    radial_velocities = [
        coord.radial_velocity.value for coord in cat_prop["new_coords"]
    ]

    # Create combined SkyCoord
    out_cat["coords"] = SkyCoord(
        ra=ras * u.deg,
        dec=decs * u.deg,
        distance=distances * u.pc,
        pm_ra_cosdec=pm_ra_cosdec * u.mas / u.yr,
        pm_dec=pm_dec * u.mas / u.yr,
        radial_velocity=radial_velocities * u.km / u.s,
        frame="icrs",
    )

    out_cat["source_id"] = np.asarray([f"Gaia DR3 {i}" for i in cat_prop["source_id"]])

    # Apply limit if necessary
    if limit is not None:
        out_cat = {k: v[:limit] for k, v in out_cat.items()}

    return out_cat


@lru_cache
def get_planets(
    #    coord: SkyCoord,
    ra: Union[float, None] = None,
    dec: Union[float, None] = None,
    name: Union[str, None] = None,
    radius: u.Quantity = 20 * u.arcsecond,
    attrs: List = [],
    # attrs: List = ["pl_orbper", "pl_tranmid", "pl_trandur", "pl_trandep"],
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
        if name is not None:
            planets_tab = NasaExoplanetArchive.query_object(name, table="pscomppars")
        elif ra is not None and dec is not None:
            planets_tab = NasaExoplanetArchive.query_region(
                table="pscomppars",
                coordinates=SkyCoord(ra, dec, unit=u.deg),
                radius=radius,
            )
        else:
            raise ValueError
        if len(planets_tab) != 0:
            # if len(attrs) == 0:
            #     attrs = planets_tab.keys()
            # else:
            #     attrs = List[attrs]
            # planets = {
            #     letter: {
            #         attr: planets_tab[planets_tab["pl_letter"] == letter][attr][
            #             0
            #         ]  # .unmasked
            #         for attr in attrs
            #     }
            #     for letter in planets_tab["pl_letter"]
            # }
            # planets = planets_tab.to_pandas()
            # planets = planets.to_dict(orient='records')
            # print(planets)

            # There's an error in the NASA exoplanet archive units that makes duration "days" instead of "hours"
            # for planet in planets:
            #     if "pl_trandur" in planets[planet].keys():
            #         planets[planet]["pl_trandur"] = (
            #             planets[planet]["pl_trandur"].value * u.hour
            #         )
            if planets_tab["pl_trandur"].unit == u.day:
                planets_tab["pl_trandur"] = planets_tab["pl_trandur"].value * u.hour
                planets_tab["pl_trandurerr1"] = (
                    planets_tab["pl_trandurerr1"].value * u.hour
                )
                planets_tab["pl_trandurerr2"] = (
                    planets_tab["pl_trandurerr2"].value * u.hour
                )

            if len(attrs) != 0:
                planets_tab = planets_tab[attrs]

        else:
            planets_tab = QTable()

    return planets_tab


@lru_cache
def get_citation(bibcode):
    """Goes to NASA ADS and webscrapes the bibtex citation for a given bibcode"""
    # d = requests.get(f"https://ui.adsabs.harvard.edu/abs/{bibcode}/exportcitation")
    # soup = BeautifulSoup(d.content, "html.parser")
    # return soup.find("textarea").text
    bibcode = bibcode.replace("%26", "&")

    TOKEN = "hAEaZKkVUdJ2sKXHKVUM7WfhyOwbq7uls1p4poUX"
    payload = {"bibcode": [bibcode], "sort": "first_author asc"}
    results = requests.post(
        "https://api.adsabs.harvard.edu/v1/export/bibtex",
        headers={"Authorization": "Bearer " + TOKEN},
        data=json.dumps(payload),
    )
    return results.json()["export"]
