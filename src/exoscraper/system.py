"""Classes for working with Targets"""

from typing import Union

import astropy.units as u
import numpy as np

# import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.time import Time
import requests

from .query import get_planets, get_SED, get_sky_catalog, get_offline_star_catalog
from .planet import Planet
from .star import Star
from .utils import get_batman_model, handle_gaiaoffline_files


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
        name: Union[str, None] = None,
        ra: Union[u.Quantity, None] = None,
        dec: Union[u.Quantity, None] = None,
        coord: Union[SkyCoord, None] = None,
        logg: Union[u.Quantity, None] = None,
        teff: Union[u.Quantity, None] = None,
        bmag: Union[u.Quantity, None] = None,
        jmag: Union[u.Quantity, None] = None,
        time: Time = Time.now(),
        offline: bool = True,
    ):
        """Ensures quantity conventions, generates Planet and Star classes, and validates input"""
        if all(x is None for x in [name, ra, dec, coord]):
            raise ValueError("Coordinate or name must be provided!")
        self.name, self.coord, self.bmag, self.jmag = (name, coord, bmag, jmag)
        # self.ra, self.dec = u.Quantity(ra, u.deg), u.Quantity(dec, u.deg)
        # self.teff = u.Quantity(teff, u.K)
        # self.logg = u.Quantity(logg)
        self.ra = ra
        self.dec = dec
        self.teff = teff
        self.logg = logg
        self.coord = coord
        self.offline = offline

        # Processing RA and Dec input
        if self.ra is None and self.dec is None:
            if coord is None:
                self.coord = SkyCoord.from_name(name)
            self.ra, self.dec = u.Quantity(self.coord.ra, u.deg), u.Quantity(
                self.coord.dec, u.deg
            )

        # Fetching Gaia DR3 values
        if self.offline:
            self.sky_cat = get_offline_star_catalog(
                self.ra, self.dec, limit=1, time=time
            )
        else:
            self.sky_cat = get_sky_catalog(self.ra, self.dec, limit=1, time=time)
        self.coord = self.sky_cat["coords"]

        # Fetching any planets from the system
        if self.name is not None:
            sys_info = get_planets(name=self.name)
        else:
            sys_info = get_planets(self.ra.value, self.dec.value)

        if len(sys_info.columns) == 0:
            if "source_id" in self.sky_cat.keys():
                self.sky_cat["hostname"] = self.sky_cat.pop("source_id")
            self.sys_info = QTable(self.sky_cat)
            self.sys_info["sy_snum"] = [1]
        else:
            self.sys_info = sys_info
            self.sys_info["ra"] = self.ra
            self.sys_info["dec"] = self.dec
            self.sys_info["sy_pmra"] = self.coord.pm_ra_cosdec
            self.sys_info["sy_pmdec"] = self.coord.pm_dec
            self.sys_info["coord_epoch"] = time

        # Loop through the unique hostnames in the query and make Star objects out of them
        # Right now everything is funneled through an Exo Archive query which does not capture
        # multi-star systems and treats each individual star in the system as a single star.
        # Maybe change the star query to something else and query Exo Archive within Star?
        self.stars = []
        for host in np.unique(self.sys_info["hostname"]):
            self.stars.append(
                Star(self.sys_info[self.sys_info["hostname"] == str(host)])
            )

        # Will need to expand this to include binaries eventually
        if self.sys_info[0]["sy_snum"] == 1:
            self.__dict__.update(Star(self.sys_info).__dict__)

            if len(sys_info.columns) > 0:
                for i in range(len(self.sys_info)):
                    setattr(
                        self,
                        str(self.sys_info["pl_letter"][i]),
                        Planet(self.sys_info[i]),
                    )

        # Loop through hostnames in query and assign their variables letter names
        st_letters = ["A", "B", "C", "D", "E", "F"]
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
    def from_gaia(coord: Union[str, SkyCoord], time=Time.now(), offline=True):
        name = None
        if isinstance(coord, str):
            name = coord
            coord = SkyCoord.from_name(coord)
        elif not isinstance(coord, SkyCoord):
            raise ValueError("`coord` must be a `SkyCoord` or a name string.")
        if offline:
            handle_gaiaoffline_files()
            cat = get_offline_star_catalog(coord.ra, coord.dec, time=time)
        else:
            try:
                cat = get_sky_catalog(
                    coord.ra, coord.dec, radius=5 * u.arcsecond, limit=1, time=time
                )
            except TimeoutError:
                print("TimeoutError: Trying gaiaoffline query")
                handle_gaiaoffline_files()
                cat = get_offline_star_catalog(coord.ra, coord.dec, time=time)
            except requests.exceptions.HTTPError as http_err:
                print(f"HTTP error occurred: {http_err}. Trying offline query.")
                handle_gaiaoffline_files()
                cat = get_offline_star_catalog(coord.ra, coord.dec, time=time)
        if name is None:
            name = cat["source_id"][0]
        return System(
            name=name,
            ra=cat["coords"][0].ra,
            dec=cat["coords"][0].dec,
            logg=cat["logg"][0],
            teff=cat["teff"][0],
            bmag=cat["bmag"][0],
            jmag=cat["jmag"][0],
            coord=cat["coords"][0],
            time=time,
        )

    @staticmethod
    def from_TIC(coord: Union[str, SkyCoord]):
        raise NotImplementedError

    @staticmethod
    def from_name(name: str, time=Time.now()):
        if not isinstance(name, str):
            raise ValueError("`name` must be a `string`.")
        else:
            return System.from_gaia(coord=name, time=time)

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
    def transit_model(self) -> np.ndarray:
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
        t0: list = [],
        period: list = [],
        ror: list = [],
        dor: list = [],
        inc: list = [],
        ecc: list = [],
        periastron: list = [],
        limb_dark: str = "uniform",
        u: list = [],
        iterations: int = 1,
        seed: Union[int, None] = None,
        median_flag: bool = False,
        vals_out: bool = False,
        **kwargs,
    ):
        """
        Passes known system information to exoplanet to generate BATMAN model. Samples a single
        iteration of the time series by default. If values are not specified for each parameter,
        the values listed for each planet on the NASA Exoplanet Archive and their associated errors
        will be used.

        Parameters
        ----------
        time : np.array
            Timestamps at which to model the timeseries.
        t0 : list
            Mid-transit times to be used for each planet in the timeseries model. Unit should be the
            same as is used in `time`.
        period : list
            Periods to be used for each planet in the timeseries model. Unit should be the same as is used
            in `time`.
        ror : list
            Ratio of planetary radius to stellar radius for each planet in the time series model.
        dor : list
            Ratio of semimajor axis to the radius of the host star for each planet in the time series model.
        inc : list
            Orbital inclination in degrees for each planet in the time series model.
        ecc : list
            Orbital eccentricity for each planet in the time series model.
        periastron : list
            Longitude of periastron in degrees for each planet in the time series model. Default is 90
            degrees.
        limb_dark : str
            Limb darkening model to use in modeling the planet transits in the time series. Default is
            "uniform".
        u : list
            Limb darkening coefficients corresponding to the limb darkening model specified in `limb_dark`.
        iterations : int
            Number of timeseries models to generate for the system. This will determine how many parameter
            draws are performed if `median_flag` is False and any parameters are left unspecified.
        seed : int
            Seed value for the draws from the parameter distributions if `median_flag` is False.
        median_flag : bool
            Flag determining if median values for each parameter are used. If set to True, the median values
            for each parameter will be used. If set to False, values will be drawn from the distributions
            of each parameter.
        vals_out : bool
            Flag determining whether the parameter values used in each timeseries model will be returned. If
            set to True, the function will return two outputs.
        """
        par_names = [
            "pl_tranmid",
            "pl_orbper",
            "pl_ratror",
            "pl_ratdor",
            "pl_orbincl",
            "pl_orbeccen",
            "pl_orblper",
        ]
        pars_in = {
            "t0": t0,
            "period": period,
            "ror": ror,
            "dor": dor,
            "inc": inc,
            "ecc": ecc,
            "peri": periastron,
        }
        vars = {
            "t0": np.ones((iterations, len(self.planets))),
            "period": np.ones((iterations, len(self.planets))),
            "ror": np.ones((iterations, len(self.planets))),
            "dor": np.ones((iterations, len(self.planets))),
            "inc": np.ones((iterations, len(self.planets))),
            "ecc": np.ones((iterations, len(self.planets))),
            "peri": np.ones((iterations, len(self.planets))),
        }
        timeseries = np.zeros((iterations, len(time)))

        if median_flag:
            for p, par in enumerate(vars.keys()):
                if len(pars_in[par]) == len(self.planets):
                    vars[par] *= pars_in[par]
                else:
                    vars[par] *= [
                        getattr(self[0][t], par_names[p]).value
                        for t in range(len(self.planets))
                    ]
        else:
            vars["peri"] = [
                getattr(self[0][t], "pl_orblper").value
                for t in range(len(self.planets))
            ]
            for p, par in enumerate(vars.keys()):
                if len(pars_in[par]) == 0:
                    vars[par] *= np.array(
                        [
                            getattr(self[0][t], par_names[p])
                            .distribution.sample(seed=seed, size=iterations)
                            .value
                            for t in range(len(self.planets))
                        ]
                    ).T
                elif len(pars_in[par]) == len(self.planets):
                    vars[par] *= pars_in[par]
                else:
                    raise ValueError(
                        "Number of parameter values provided must match number of planets"
                    )

        for i in range(iterations):
            flux = np.zeros(len(time))

            for n, pl in enumerate(self.planets):
                model, params = get_batman_model(
                    time=time,
                    t0=vars["t0"][i][n],
                    per=vars["period"][i][n],
                    ror=vars["ror"][i][n],
                    dor=vars["dor"][i][n],
                    inc=vars["inc"][i][n],
                    ecc=vars["ecc"][i][n],
                    periastron=vars["peri"][i][n],
                    limb_dark=limb_dark,
                    u=u,
                    params_out=True,
                    **kwargs,
                )
                setattr(pl, "model", model)
                flux += model.light_curve(params)

            timeseries[i] += flux - len(self.planets) + 1

        if vals_out:
            return timeseries, vars
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
