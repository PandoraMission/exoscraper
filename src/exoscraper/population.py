"""Classes for working with Populations"""

from pandas import DataFrame


class Population(DataFrame):
    """Work with the exoplanet archive population?"""

    def __init__(self):
        # names: list
        # ra: array of u.Quantity
        # dec: array of u.Quantity
        # boundaries: dict

        # get_params(ra, dec, names, boundaries)
        raise NotImplementedError

    # Go get the data

    # Some sort of plotting tools or diagnostics

    @property
    def to_TargetSet(self):
        # converts the current population selection into a TargetSet
        raise NotImplementedError
