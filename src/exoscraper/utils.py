"""File to hold a bunch of utility functions"""

import numpy as np
import batman
import pandas as pd

from gaiaoffline import Gaia, populate_gaiadr3, populate_tmass, populate_tmass_xmatch

from .query import get_citation


def get_ref_dict(tab):
    """Parses the NExSci table for a list of references"""
    cols = [c for c in tab.columns if "reflink" in c]
    refs = np.unique(tab[cols])[0]
    result = {
        ref.split(">")[1]
        .split("</a")[0]
        .strip(): ref.split("href=")[1]
        .split(" target=ref")[0]
        for ref in refs
        if ref != ""
    }
    for key, item in result.items():
        if "ui.adsabs" in item.lower():
            result[key] = get_citation(item.split("abs/")[1].split("/")[0])
    return result


def get_batman_model(
    time: np.array,
    t0: float,
    per: float,
    ror: float,
    dor: float,
    inc: float = 90.0,
    ecc: float = 0.0,
    periastron: float = 90.0,
    limb_dark: str = "uniform",
    u: list = [],
    params_out: bool = False,
    **kwargs,
):
    """Generates a batman model of the exoplanet orbit"""
    params = batman.TransitParams()
    params.t0 = t0
    params.per = per  # days
    params.rp = ror  # stellar radius
    params.a = dor  # stellar radius
    params.inc = inc  # degrees
    params.ecc = ecc
    params.w = periastron  # longitude of periastron
    params.limb_dark = limb_dark  # limb darkening model
    params.u = u  # limb darkening parameters

    model = batman.TransitModel(params, time, **kwargs)

    if params_out:
        return model, params
    else:
        return model


def handle_gaiaoffline_files():
    """
    Function to check the local `gaiaoffline` install and add the relevant files if necessary.
    """
    tracker_table_names = Gaia().file_tracker_table_names
    track_fracs = {}
    for tracker_table_name in tracker_table_names:
        df = pd.read_sql_query(f"""SELECT * FROM {tracker_table_name}""", Gaia().conn)
        frac = df.status.isin(["completed"]).sum() / len(df)
        track_fracs.update({str(tracker_table_name): frac})

    if ("file_tracking_gaiadr3" not in track_fracs.keys()) or track_fracs[
        "file_tracking_gaiadr3"
    ] < 1.0:
        print("Warning: This may take a while.")
        print("Populating Gaia DR3 files for gaiaoffline")
        populate_gaiadr3()

    if ("file_tracking_tmass_xmatch" not in track_fracs.keys()) or track_fracs[
        "file_tracking_tmass_xmatch"
    ] < 1.0:
        print("Warning: This may take a while.")
        print("Populating 2MASS Xmatch files for gaiaoffline")
        populate_tmass_xmatch()

    if ("file_tracking_tmass" not in track_fracs.keys()) or track_fracs[
        "file_tracking_tmass"
    ] < 1.0:
        print("Warning: This may take a while.")
        print("Populating 2MASS files for gaiaoffline")
        populate_tmass()


def calc_separation(ref_ra, ref_dec, ra, dec):
    """
    Function to calculate the separation between two RA/Dec coordinate sets.

    ref_ra: float
        RA of reference coordinate set.
    ref_dec: float
        Dec of reference coordinate set.
    ra: float
        RA of target coordinate set.
    dec: float
        Dec of target coordinate set.
    """
    ra1 = np.radians(ref_ra).astype(np.float64)
    ra2 = np.radians(ra).astype(np.float64)
    dec1 = np.radians(ref_dec).astype(np.float64)
    dec2 = np.radians(dec).astype(np.float64)

    numerator = (
        np.sin((dec2 - dec1) / 2) ** 2
        + np.cos(dec1) * np.cos(dec2) * np.sin((ra2 - ra1) / 2) ** 2
    )
    sep = 2 * np.arctan2(numerator**0.5, (1 - numerator) ** 0.5)

    return sep
