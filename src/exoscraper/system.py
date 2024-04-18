"""Classes for working with Targets"""
from typing import Union

import astropy.units as u
import numpy as np
# import pandas as pd
import lightkurve as lk
from astropy.coordinates import SkyCoord

from .query import get_planets, get_SED, get_sky_catalog
from .planet import Planet
from .star import Star
from .utils import get_batman_model


class System(object):
    # name: str = None
    # ra: u.Quantity = None
    # dec: u.Quantity = None
    # logg: u.Quantity
    # teff: u.Quantity
    # bmag: u.Quantity
    # jmag: u.Quantity
    # coord: SkyCoord = None
    """DOCSTRING"""

    def __init__(
            self,
            name: str | None = None,
            ra: u.Quantity | None = None,
            dec: u.Quantity | None = None,
            coord: SkyCoord | None = None,
            # logg: u.Quantity | None = None,
            # teff: u.Quantity | None = None,
            bmag: u.Quantity | None = None,
            jmag: u.Quantity | None = None,
            ):
        """Ensures quantity conventions, generates Planet and Star classes, and validates input"""
        if all(x is None for x in [name, ra, dec, coord]):
            raise ValueError
        self.name, self.coord, self.bmag, self.jmag = (name, coord, bmag, jmag)
        self.ra, self.dec = u.Quantity(ra, u.deg), u.Quantity(dec, u.deg)
        # self.teff = u.Quantity(teff, u.K)
        # self.logg = u.Quantity(logg)
        if self.name is not None:
            self.sys_info = get_planets(name=self.name)
        elif self.ra is not None and self.dec is not None:
            self.sys_info = get_planets(self.ra.value, self.dec.value)

        # Loop through the unique hostnames in the query and make Star objects out of them
        # Right now everything is funneled through an Exo Archive query which does not capture
        # multi-star systems and treats each individual star in the system as a single star.
        # Maybe change the star query to something else and query Exo Archive within Star?
        self.stars = []
        for host in np.unique(self.sys_info['hostname']):
            self.stars.append(Star(self.sys_info[self.sys_info['hostname'] == str(host)]))

        # Will need to expand this to include binaries eventually
        if self.sys_info[0]['sy_snum'] == 1:
            self.__dict__.update(Star(self.sys_info).__dict__)

            for i in range(len(self.sys_info)):
                setattr(self, str(self.sys_info['pl_letter'][i]), Planet(self.sys_info[i]))

        # Loop through hostnames in query and assign their variables letter names
        st_letters = ['A', 'B', 'C', 'D', 'E', 'F']
        for i in range(len(self.stars)):
            setattr(self, st_letters[i], self.stars[0])

        # for i in range(len(self.sys_info)):
        #     # setattr(self, 'planet' + str(self.sys_info['pl_letter'][i]), Planet(self.sys_info[i]))
        #     setattr(self, str(self.sys_info['pl_letter'][i]), Planet(self.sys_info[i]))

        return

    def __repr__(self):
        return f"{self.name} [{self.ra}, {self.dec}]"

    def _repr_html_(self):
        return f"{self.name} ({self.ra._repr_latex_()},  {self.dec._repr_latex_()})"

    def __getitem__(self, index):
        return self.stars[index]

    @staticmethod
    def from_gaia(coord: Union[str, SkyCoord]):
        name = None
        if isinstance(coord, str):
            name = coord
            coord = SkyCoord.from_name(coord)
        elif not isinstance(coord, SkyCoord):
            raise ValueError("`coord` must be a `SkyCoord` or a name string.")
        cat = get_sky_catalog(coord.ra, coord.dec, radius=5 * u.arcsecond, limit=1)
        if name is None:
            name = cat["source_id"][0]
        return System(
            name=name,
            ra=cat["coords"][0].ra,
            dec=cat["coords"][0].dec,
            # logg=cat["logg"][0],
            # teff=cat["teff"][0],
            bmag=cat["bmag"][0],
            jmag=cat["jmag"][0],
            coord=cat["coords"][0],
        )

    @staticmethod
    def from_TIC(coord: Union[str, SkyCoord]):
        raise NotImplementedError

    @staticmethod
    def from_name(name: str):
        if not isinstance(name, str):
            raise ValueError("`name` must be a `string`.")
        else:
            return System.from_gaia(coord=name)

    @property
    def SED(self) -> dict:
        """Returns a dictionary containing the SED of the target from Vizier

        Uses a default radius of 3 arcseconds.
        """
        return get_SED((self.ra.value, self.dec.value), radius=3 * u.arcsecond)

    @property
    def lightcurve(self):
        # go get the TESS data something like
        # return get_timeseries(self.ra, self.dec)
        raise NotImplementedError

    @property
    def bibliography(self):
        # go get references from NASA ADS
        # return get_bibliography(self.name)
        raise NotImplementedError

    @property
    def transit_times(self, time: u.Quantity):
        # calculate future transit times
        # return array of future transit times out to a specified time
        raise NotImplementedError

    @property
    def transit_model(self) -> lk.LightCurve:
        # generate transit model using lightkurve
        # return LightCurve object
        raise NotImplementedError

    @property
    def noise_model(self, mission: str) -> float:
        # generate noise model for target for specified mission
        # this could be wrapped into transit model
        # this also may not be worth including since I selfishly need it for
        #  pandora-target and so added it here lol
        # returns LightCurve object
        raise NotImplementedError

    def sample_timeseries(
            self,
            time,
            t0: list | None = None,
            period: list | None = None,
            ror: list | None = None,
            dor: list | None = None,
            inc: list | None = None,
            ecc: list | None = None,
            periastron=90,
            limb_dark: str = "uniform",
            u: list = [],
            iterations: int = 1,
            seed: int | None = None,
            median_flag: bool = False,
            vars_out: bool = False,
            **kwargs,
    ):
        """
        Passes known system information to exoplanet to generate BATMAN model. Samples a single
        iteration of the time series by default.
        """
        pars = ['pl_tranmid', 'pl_orbper', 'pl_ratror', 'pl_ratdor', 'pl_orbincl', 'pl_orbeccen']
        timeseries = np.zeros((iterations, len(time)))
        vars_list = np.zeros((iterations, len(pars)))
        for i in range(iterations):
            # build input arrays for all variables that are not user-provided
            for m, var in enumerate([t0, period, ror, dor, inc, ecc]):
                if median_flag:
                    var = [getattr(self[0][t], pars[m]).value for t in range(len(self.planets))]
                if var is None:
                    var = [getattr(self[0][t], pars[m]).distribution.sample(seed=seed).value for t in range(len(self.planets))]
                if vars_out:
                    vars_list[i][m] = var

            flux = np.zeros(len(time))

            for n, pl in enumerate(self.planets):
                model, params = get_batman_model(
                    time=time,
                    t0=t0[n][i],
                    per=period[n][i],
                    ror=ror[n][i],
                    dor=dor[n][i],
                    inc=inc[n][i],
                    ecc=ecc[n][i],
                    periastron=periastron,
                    limb_dark=limb_dark,
                    u=u,
                    params_out=True,
                    **kwargs
                )
                setattr(pl, "model", model)
                flux += model.light_curve(params)

            timeseries[i] += flux

        if vars_out:
            return timeseries, vars_list
        else:
            return timeseries

    def motto(self):
        # mottos = ["If it's online, we can scrape it!",
        #           "Leave the scraping to us!",
        #           "Scraping the internet clean!",
        #           "Scrape it or leave it!",
        #           "Consider it scraped!"]

        # print(random.choice(mottos))
        raise NotImplementedError


# class SystemSet(object):
#     """A class to hold many Target classes"""

#     def __init__(self, targets):
#         self.targets = targets

#     def __iter__(self):
#         raise NotImplementedError

#     def __len__(self):
#         raise NotImplementedError

#     def __repr__(self):
#         raise NotImplementedError

#     @staticmethod
#     def from_names(coords: Union[str, SkyCoord]):
#         raise NotImplementedError

#     def to_csv(self, output: str):
#         """Produces csv file with all the targets in the TargetSet and saves to output
#         Parameters
#         ----------
#         output: string
#             csv output location and desired file name

#         Return
#         ------
#         csv file
#             file containing list of planets and select parameters for all targets in TargetSet
#         """

#         # Initialize DataFrame to fill with Targets from TargetSet
#         targets_df = pd.DataFrame(
#             [],
#             columns=[
#                 "Planet Name",
#                 "Star Name",
#                 "Star SkyCoord",
#                 "Planet Transit Epoch (BJD_TDB-2400000.5)",
#                 "Planet Transit Epoch Uncertainty",
#                 "Period (day)",
#                 "Period Uncertainty"
#                 "Transit Duration (hrs)",
#             ],
#         )

#         # Pull data from TargetSet to fill DataFrame
#         # for target in targets:

#         # Save DataFrame to csv
#         targets_df = targets_df.sort_values(by=["Planet Name"]).reset_index(drop=True)
#         targets_df.to_csv((output), sep=",", index=False)
#         raise NotImplementedError
