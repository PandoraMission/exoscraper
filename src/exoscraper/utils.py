"""File to hold a bunch of utility functions"""

import numpy as np
import batman

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
