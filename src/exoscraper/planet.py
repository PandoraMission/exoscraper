"""Classes to work specifically with planets"""

from dataclasses import dataclass

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.constants import G
from astropy.table import QTable

from .query import get_citation
from .utils import get_ref_dict
from .distribution import NormalDistribution, UniformDistribution


@dataclass
class Planet(object):
    """Helper class to hold planet information from NExSci. This class only holds and prints information,
    it doesn't calculate anything.
    """

    def __init__(
        self,
        params: QTable,
    ):
        self.hostname = params["hostname"]
        self.letter = params["pl_letter"]
        self._tab = params
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
        # Error on the archive that the unit is days when it should be hours!
        if self.pl_trandur.unit == u.day:
            # self.pl_trandur /= 24  # I don't think the value needs to be changed, just the unit
            self.pl_trandur = self.pl_trandur.value * u.hour

        for c in self._tab.columns:
            if c.endswith("reflink"):
                attr = getattr(self, "_".join(c.split("_")[:-1]))
                if isinstance(attr, u.Quantity):
                    if self._tab[c] != "":
                        ref = str(
                            self._tab[c].split("href=")[1].split(" target=ref")[0]
                        )
                        if "ui.adsabs" in ref.lower():
                            # print(ref.split("abs/")[1].split("/")[0])
                            # print(ref)
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
                        elif "calculated value" in ref.lower():
                            setattr(attr, "reference", "Archive Calculation")
                        elif "exofop" in ref.lower():
                            setattr(attr, "reference", "ExoFOP")
                        elif "arxiv.org" in ref.lower():
                            ref = self._tab[c].split("href=")[1].split(" target=ref")[0]
                            setattr(attr, "reference", ref)
                            setattr(
                                attr,
                                "reference_name",
                                self._tab[c].split("=ref>")[1].split("</a")[0],
                            )
                            setattr(attr, "reference_link", ref)
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
        self._make_distributions()
        self.references = get_ref_dict(self._tab)
        self.acknowledgements = [
            "This research has made use of the NASA Exoplanet Archive, which is operated by the"
            "  California Institute of Technology, under contract with the National Aeronautics "
            "and Space Administration under the Exoplanet Exploration Program."
        ]
        return

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
        elif self.pl_orbsmax.unit is not u.AU:
            self.pl_orbsmax._unit = u.AU

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

    def _make_distributions(self):
        """Hidden function to build distributions for each parameter"""
        for c in self._tab.columns:
            if c.endswith("err1"):
                attr = getattr(self, c[:-4])
                if isinstance(attr, u.Quantity):
                    if self._tab[c] != "" and not np.isnan(attr.value):
                        reflink = None
                        if ((c[:-4] + "_reflink") in self._tab.columns) & (
                            self._tab[c[:-4] + "_reflink"] != ""
                        ):
                            reflink = self._tab[c[:-4] + "_reflink"]
                        else:
                            setattr(attr, "reference", None)
                        if hasattr(attr, "err1"):
                            err = max(abs(attr.err1.value), abs(attr.err2.value))
                        else:
                            err = attr.err
                        setattr(
                            attr,
                            "distribution",
                            NormalDistribution(
                                attr.value,
                                err,
                                name=c[:-4],
                                unit=str(attr.unit),
                                reference=reflink,
                            ),
                        )
            if c.endswith("lim"):
                if self._tab[c] == 1:
                    attr = getattr(self, c[:-3])
                    reflink = None
                    if (c[:-3] + "_reflink") in self._tab.columns:
                        reflink = attr.reference_link
                    setattr(
                        attr,
                        "distribution",
                        UniformDistribution(
                            0,
                            self._tab[c[:-3]],
                            name=c[:-3],
                            unit=str(attr.unit),
                            reference=reflink,
                        ),
                    )

    def __repr__(self):
        return self.hostname + self.letter

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


# class Planets(object):
#     """Special class to hold many planets in one system"""

#     def __init__(self, hostname: str):
#         self.hostname = hostname
#         tab = get_nexsci_tab(hostname)
#         if len(tab) == 0:
#             raise ValueError("No planets found")
#         self.letters = np.unique(list(tab["pl_letter"]))
#         self.planets = [Planet(self.hostname, letter) for letter in self.letters]
#         self._tab = tab
#         self._cols = [
#             c
#             for c in list(self._tab.columns)
#             if not (
#                 c.endswith("err1")
#                 | c.endswith("err2")
#                 | c.endswith("reflink")
#                 | c.endswith("lim")
#                 | c.endswith("str")
#             )
#         ]
#         _ = [
#             setattr(self, attr, [getattr(planet, attr) for planet in self])
#             for attr in self._cols
#             if attr.startswith("pl")
#         ]
#         _ = [
#             setattr(self, attr, getattr(self[0], attr))
#             for attr in self._cols
#             if attr.startswith("st")
#         ]
#         _ = [
#             setattr(self, attr, getattr(self[0], attr))
#             for attr in self._cols
#             if attr.startswith("sy")
#         ]
#         self.acknowledgements = []
#         self.references = self[0].references
#         for planet in self:
#             self.references.update(planet.references)
#             [self.acknowledgements.append(a) for a in planet.acknowledgements]
#         self.acknowledgements = list(np.unique(self.acknowledgements))

#     def __len__(self):
#         return len(self.letters)

#     def __repr__(self):
#         return (
#             f"{self.hostname} System ({len(self)} Planet{'s' if len(self) > 1 else ''})"
#         )

#     def __getitem__(self, idx):
#         if isinstance(idx, int):
#             return self.planets[idx]
#         elif isinstance(idx, str):
#             if idx in self.letters:
#                 return self.planets[np.where(self.letters == idx)[0][0]]
#             else:
#                 raise ValueError(f"No planet `{idx}` in the {self.hostname} system.")
#         else:
#             raise ValueError(f"Can not parse `{idx}` as a planet.")

#     @property
#     def StarParametersTable(self):
#         return self[0].StarParametersTable

#     @property
#     def StarParametersTableLatex(self):
#         return self[0].StarParametersTableLatex

#     @property
#     def PlanetsParametersTable(self):
#         dfs = [
#             planet.PlanetParametersTable[["Value", "Reference"]].rename(
#                 {"Value": planet.name}, axis="columns"
#             )
#             for planet in self
#         ]
#         return pd.concat([self[0].PlanetParametersTable[["Description"]], *dfs], axis=1)

#     @property
#     def PlanetsParametersTableLatex(self):
#         print(
#             self.PlanetsParametersTable.to_latex(
#                 caption=f"Planet Parameters for {self.hostname} Planets",
#                 label="tab:planetparams",
#             )
#         )
