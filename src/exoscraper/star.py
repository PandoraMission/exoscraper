"""Defines the dataclass to hold stellar information"""

from dataclasses import dataclass
from typing import Union, Dict

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.table import QTable

from .utils import get_ref_dict
from .query import get_citation
from .planet import Planet


@dataclass
class Star(object):
    """Helper class to hold information about stars in a system. The only calculations performed
    in this class are those necessary to clean the data.
    """

    def __init__(
        self,
        params: Union[QTable, Dict],
    ):
        if type(params) is Dict:
            params = QTable(params)

        self.name = params["hostname"][0]

        if "pl_letter" in params.columns:
            self.planets = []
            for i, letter in enumerate(np.unique(params["pl_letter"])):
                self.planets.append(
                    Planet(params[params["pl_letter"] == str(letter)][0])
                )
                setattr(self, str(letter), self.planets[i])

        good_inds = [col for col in params.columns if "pl_" not in col]
        self._tab = params[0][good_inds]

        _ = [
            (
                setattr(
                    self,
                    c,
                    (
                        np.ma.MaskedArray(self._tab[c]).filled(np.nan)
                        if isinstance(self._tab[c], u.Quantity)
                        else u.Quantity(self._tab[c])
                    ),
                )
                if isinstance(self._tab[c], (u.Quantity, float, int))
                else setattr(self, c, self._tab[c])
            )
            for c in list(self._tab.columns)
            if not (
                c.endswith("err1")
                | c.endswith("err2")
                | c.endswith("reflink")
                | c.endswith("lim")
                | c.endswith("str")
            )
        ]

        for c in self._tab.columns:
            if c.endswith("reflink"):
                attr = getattr(self, "_".join(c.split("_")[:-1]))
                if isinstance(attr, u.Quantity):
                    if self._tab[c] != "":
                        ref = self._tab[c].split("href=")[1].split(" target=ref")[0]
                        if "ui.adsabs" in ref.lower():
                            ref = get_citation(ref.split("abs/")[1].split("/")[0])
                            setattr(attr, "reference", ref)
                            setattr(
                                attr, "reference_name", ref.split("{")[1].split(",")[0]
                            )
                            setattr(
                                attr,
                                "reference_link",
                                ref.split("adsurl = {")[1].split("}")[0],
                            )
            if c.endswith("err1"):
                attr = getattr(self, c[:-4])
                if isinstance(attr, u.Quantity):
                    if self._tab[c] != "":
                        setattr(attr, "err1", u.Quantity(self._tab[c[:-4] + "err1"]))
                        setattr(attr, "err2", u.Quantity(self._tab[c[:-4] + "err2"]))
                        setattr(
                            attr,
                            "err",
                            u.Quantity(
                                [
                                    self._tab[c[:-4] + "err1"],
                                    -self._tab[c[:-4] + "err2"],
                                ]
                            ).mean(),
                        )
                    else:
                        setattr(self._tab[c], "err", np.nan * self._tab[c].unit)
            if c.endswith("lim"):
                # Any "limit" parameters need to be set to nans
                if self._tab[c] == 1:
                    attr = getattr(self, c[:-3])
                    attr *= np.nan
                    for e in ["err1", "err2", "err"]:
                        if hasattr(attr, e):
                            setattr(attr, e, getattr(attr, e) * np.nan)

        if any([c.endswith("reflink") for c in self._tab.columns]):
            self.references = get_ref_dict(self._tab)
        self.acknowledgements = [
            "This research has made use of the NASA Exoplanet Archive, which is operated by the"
            "  California Institute of Technology, under contract with the National Aeronautics "
            "and Space Administration under the Exoplanet Exploration Program."
        ]
        return

    def __repr__(self):
        return self.name

    def __getitem__(self, index):
        return self.planets[index]

    @property
    def StarParametersTable(self):
        d = pd.DataFrame(columns=["Value", "Description", "Reference"])
        for key, symbol, desc in zip(
            ["st_rad", "st_mass", "st_age", "st_logg"],
            ["R", "M", "Age", "logg"],
            ["Stellar Radius", "Stellar Mass", "Stellar Age", "Stellar Gravity"],
        ):
            attr = getattr(self, key)
            if np.isfinite(attr):
                d.loc[symbol, "Value"] = "{0}^{{{1}}}_{{{2}}}".format(
                    attr.to_string(format="latex"), attr.err1.value, attr.err2.value
                )
                d.loc[symbol, "Description"] = desc
                d.loc[symbol, "Reference"] = f"\\cite{{{attr.reference_name}}}"
        return d

    @property
    def StarParametersTableLatex(self):
        print(
            self.StarParametersTable.to_latex(
                caption=f"Stellar Parameters for {self.hostname}",
                label="tab:stellarparams",
            )
        )
