"""Classes to work specifically with planets and planet host stars"""
from functools import lru_cache

import astropy.units as u
import numpy as np
import pandas as pd
import requests
from astropy.constants import G
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from bs4 import BeautifulSoup


@lru_cache
def get_nexsci_tab(name):
    """LRU Cached helper to get the composite entry for a given name"""
    return NasaExoplanetArchive.query_object(name, table="pscomppars")


@lru_cache
def get_citation(bibcode):
    """Goes to NASA ADS and webscrapes the bibtex citation for a give bibcode"""
    d = requests.get(f"https://ui.adsabs.harvard.edu/abs/{bibcode}/exportcitation")
    soup = BeautifulSoup(d.content, "html.parser")
    return soup.find("textarea").text


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


class Planet(object):
    """Helper class to get information from NExSci. This class only holds and prints information, it doesn't calculate anything."""

    def __init__(self, hostname: str, letter: str = "b"):
        self.hostname = hostname
        self.letter = letter
        tab = get_nexsci_tab(hostname)
        tab = tab[tab["pl_letter"] == letter]
        if len(tab) == 0:
            raise ValueError("No planet found")
        self._tab = tab[0]
        _ = [
            setattr(
                self,
                c,
                self._tab[c].filled(np.nan)
                if isinstance(self._tab[c], u.Quantity)
                else u.Quantity(self._tab[c]),
            )
            if isinstance(self._tab[c], (u.Quantity, float, int))
            else setattr(self, c, self._tab[c])
            for c in list(self._tab.columns)
            if not (
                c.endswith("err1")
                | c.endswith("err2")
                | c.endswith("reflink")
                | c.endswith("lim")
                | c.endswith("str")
            )
        ]
        # Error on the archive that the unit is days not hours...!
        if self.pl_trandur.unit == u.day:
            self.pl_trandur /= 24

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
        self._fix_eorj()
        self._fix_orbsmax()
        self._fix_ratdor()
        self._fix_ratror()
        self._fix_eqt()
        self.references = get_ref_dict(self._tab)
        self.acknowledgements = [
            "This research has made use of the NASA Exoplanet Archive, which is operated by the California"
            " Institute of Technology, under contract with the National Aeronautics and Space Administration "
            "under the Exoplanet Exploration Program."
        ]

    @property
    def name(self):
        return self.hostname + self.letter

    def _fix_eorj(self):
        if np.isfinite(self.pl_bmasse) ^ np.isfinite(self.pl_bmassj):
            if np.isfinite(self.pl_bmasse):
                self.pl_bmassj = self.pl_bmasse.to(u.jupiterMass)
                for e in ["err1", "err2", "err"]:
                    if hasattr(self.pl_bmasse, e):
                        setattr(
                            self.pl_bmassj,
                            e,
                            getattr(self.pl_bmasse, e).to(u.jupiterMass),
                        )
                self.pl_bmassj.reference = self.pl_bmasse.reference
            else:
                self.pl_bmasse = self.pl_bmassj.to(u.earthMass)
                for e in ["err1", "err2", "err"]:
                    if hasattr(self.pl_bmassj, e):
                        setattr(
                            self.pl_bmasse,
                            e,
                            getattr(self.pl_bmassj, e).to(u.earthMass),
                        )
                self.pl_bmasse.reference = self.pl_bmassj.reference

        if np.isfinite(self.pl_rade) ^ np.isfinite(self.pl_radj):
            if np.isfinite(self.pl_rade):
                self.pl_radj = self.pl_rade.to(u.jupiterRad)
                for e in ["err1", "err2", "err"]:
                    if hasattr(self.pl_rade, e):
                        setattr(
                            self.pl_radj, e, getattr(self.pl_rade, e).to(u.jupiterRad)
                        )
            else:
                self.pl_rade = self.pl_radj.to(u.earthRad)
                for e in ["err1", "err2", "err"]:
                    if hasattr(self.pl_radj, e):
                        setattr(
                            self.pl_rade, e, getattr(self.pl_radj, e).to(u.earthRad)
                        )

    def _fix_orbsmax(self):
        if not np.isfinite(self.pl_orbsmax):
            a = u.Quantity(
                (
                    ((G * self.st_mass) / (4 * np.pi) * self.pl_orbper**2) ** (1 / 3)
                ).to(u.AU)
            )
            a.err = (
                (self.st_mass.err / self.st_mass) ** 2
                + (self.pl_orbper.err / self.pl_orbper) ** 2
            ) ** 0.5 * a
            a.reference = "Calculated"
            self.pl_orbsmax = a

    def _fix_ratdor(self):
        if not np.isfinite(self.pl_ratdor):
            q = u.Quantity(self.pl_orbsmax.to(u.AU) / self.st_rad.to(u.AU))
            q.err = (
                (
                    (self.pl_orbsmax.err / self.pl_orbsmax) ** 2
                    + (self.st_rad.err / self.st_rad) ** 2
                )
            ) ** 0.5 * q
            q.reference = "Calculated"
            self.pl_ratdor = q

    def _fix_eqt(self):
        if not np.isfinite(self.pl_eqt):
            # Assume albedo is 1
            eqt = self.st_teff * np.sqrt(0.5 * 1 / self.pl_ratdor)
            eqt.err = (self.st_teff + self.st_teff.err) * np.sqrt(
                0.5 * 1 / (self.pl_ratdor - self.pl_ratdor.err)
            ) - eqt
            eqt.reference = "Calculated"
            self.pl_eqt = eqt

    def _fix_ratror(self):
        if not np.isfinite(self.pl_ratror):
            r = self.pl_rade.to(u.solRad) / self.st_rad.to(u.solRad)
            r.err = (
                (self.pl_rade.err / self.pl_rade) ** 2
                + (self.st_rad.err / self.st_rad) ** 2
            ) ** 0.5 * r
            self.pl_ratror = r

    def __repr__(self):
        return self.hostname + self.letter

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

    @property
    def PlanetParametersTable(self):
        d = pd.DataFrame(columns=["Value", "Description", "Reference"])
        for key, symbol, desc in zip(
            ["pl_radj", "pl_bmassj", "pl_orbper", "pl_tranmid"],
            ["R", "M", "P", "T_0"],
            [
                "Planet Radius",
                "Planet Mass",
                "Planet Orbital Period",
                "Planet Transit Midpoint",
            ],
        ):
            attr = getattr(self, key)
            if np.isfinite(attr):
                d.loc[symbol, "Value"] = "{0}^{{{1}}}_{{{2}}}{3}".format(
                    attr.value,
                    attr.err1.value,
                    attr.err2.value,
                    attr.unit.to_string("latex"),
                )
                d.loc[symbol, "Description"] = desc
                d.loc[symbol, "Reference"] = (
                    f"\\cite{{{attr.reference_name}}}"
                    if hasattr(attr, "reference_name")
                    else ""
                )
        return d

    @property
    def PlanetParametersTableLatex(self):
        print(
            self.PlanetParametersTable.to_latex(
                caption=f"Planet Parameters for {self.hostname + self.letter}",
                label="tab:planetparams",
            )
        )


class Planets(object):
    """Special class to hold many planets in one system"""

    def __init__(self, hostname: str):
        self.hostname = hostname
        tab = get_nexsci_tab(hostname)
        if len(tab) == 0:
            raise ValueError("No planets found")
        self.letters = np.unique(list(tab["pl_letter"]))
        self.planets = [Planet(self.hostname, letter) for letter in self.letters]
        self._tab = tab
        self._cols = [
            c
            for c in list(self._tab.columns)
            if not (
                c.endswith("err1")
                | c.endswith("err2")
                | c.endswith("reflink")
                | c.endswith("lim")
                | c.endswith("str")
            )
        ]
        _ = [
            setattr(self, attr, [getattr(planet, attr) for planet in self])
            for attr in self._cols
            if attr.startswith("pl")
        ]
        _ = [
            setattr(self, attr, getattr(self[0], attr))
            for attr in self._cols
            if attr.startswith("st")
        ]
        _ = [
            setattr(self, attr, getattr(self[0], attr))
            for attr in self._cols
            if attr.startswith("sy")
        ]
        self.acknowledgements = []
        self.references = self[0].references
        for planet in self:
            self.references.update(planet.references)
            [self.acknowledgements.append(a) for a in planet.acknowledgements]
        self.acknowledgements = list(np.unique(self.acknowledgements))

    def __len__(self):
        return len(self.letters)

    def __repr__(self):
        return (
            f"{self.hostname} System ({len(self)} Planet{'s' if len(self) > 1 else ''})"
        )

    def __getitem__(self, idx):
        if isinstance(idx, int):
            return self.planets[idx]
        elif isinstance(idx, str):
            if idx in self.letters:
                return self.planets[np.where(self.letters == idx)[0][0]]
            else:
                raise ValueError(f"No planet `{idx}` in the {self.hostname} system.")
        else:
            raise ValueError(f"Can not parse `{idx}` as a planet.")

    @property
    def StarParametersTable(self):
        return self[0].StarParametersTable

    @property
    def StarParametersTableLatex(self):
        return self[0].StarParametersTableLatex

    @property
    def PlanetsParametersTable(self):
        dfs = [
            planet.PlanetParametersTable[["Value", "Reference"]].rename(
                {"Value": planet.name}, axis="columns"
            )
            for planet in self
        ]
        return pd.concat([self[0].PlanetParametersTable[["Description"]], *dfs], axis=1)

    @property
    def PlanetsParametersTableLatex(self):
        print(
            self.PlanetsParametersTable.to_latex(
                caption=f"Planet Parameters for {self.hostname} Planets",
                label="tab:planetparams",
            )
        )
