# import pymc3 as pm
import numpy as np
import astropy.units as u
from typing import Optional


class Distribution(object):
    """Base class to work with distributions that have errors, units and references"""

    def __init__(
        self, *args, name: str, unit: str = "", reference: Optional[str] = None
    ):
        self.name = name
        self.unit = u.Unit(unit)
        self.reference = reference
        if len(args) == 0:
            raise ValueError("Can not initialize, pass in at least one argument")
        updated_args = [np.atleast_1d(args[0])]
        lens = [len(np.atleast_1d(arg)) for arg in args]
        for length, arg in zip(lens[1:], args[1:]):
            if length == 1:
                updated_args.append(np.ones_like(updated_args[0]) * (arg))
            elif length != lens[0]:
                raise ValueError("Args must be length 1 or sames arg1")
            else:
                updated_args.append(np.atleast_1d(arg))
        self.args = updated_args

    def __len__(self):
        return len(self.args[0])

    def to_string(self, **kwargs):
        return (self.args[0] * self.unit).to_string(**kwargs)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return f"{self.to_string()}"

    def _repr_latex_(self):
        return f"{self.name} {self.to_string(format='latex', subfmt='inline')}"

    def to(self, other):
        a = (self.unit).to(other)
        return type(self)(
            self.name,
            *[arg * a for arg in self.args],
            unit=other,
            reference=self.reference,
        )


class NormalDistribution(Distribution):
    @property
    def disttype(self):
        return "Normal"

    @property
    def mean(self):
        return self.args[0]

    @property
    def sigma(self):
        if len(self.args) <= 1:
            raise ValueError("Can not initialize, pass in sigma value")
        return self.args[1]

    @property
    def distribution(self):
        # return pm.Normal(self.name, self.mean, self.sigma, shape=len(self))
        raise NotImplementedError

    def sample(self, size=1, seed=None):
        return (
            np.random.default_rng(seed).normal(
                self.mean,
                self.sigma,
                size=(
                    size,
                    # len(self),
                ),
            )
            * self.unit
        )

    def to_string(self, **kwargs):
        return (
            (self.mean[0] * self.unit).to_string(**kwargs)
            if len(self) == 1
            else (self.mean * self.unit).to_string(**kwargs)
        )

    def _repr_latex_(self):
        return (
            f"{self.name} ({self.disttype}) {self.to_string(format='latex', subfmt='inline')} $\pm$ "
            + f"{(self.sigma[0] * self.unit if len(self) == 1 else self.sigma * self.unit).to_string(format='latex', subfmt='inline')}"
        )


class LogNormalDistribution(NormalDistribution):
    @property
    def disttype(self):
        return "LogNormal"

    def sample(self, size=1):
        raise NotImplementedError

    @property
    def distribution(self):
        # return pm.LogNormal(self.name, self.mean, self.sigma, shape=len(self))
        raise NotImplementedError


class UniformDistribution(Distribution):
    @property
    def disttype(self):
        return "Uniform"

    @property
    def lower(self):
        return self.args[0]

    @property
    def upper(self):
        if len(self.args) <= 1:
            raise ValueError("Can not initialize, pass in upper value")
        return self.args[1]

    @property
    def mean(self):
        return np.mean([self.lower, self.upper], axis=0)

    @property
    def diff(self):
        return (self.upper - self.lower) / 2

    @property
    def distribution(self):
        # return pm.Uniform(self.name, self.lower, self.upper, shape=len(self))
        raise NotImplementedError

    def sample(self, size=1, seed=None):
        return (
            np.random.default_rng(seed).uniform(
                self.lower,
                self.upper,
                size=(
                    size,
                    # len(self),
                ),
            )
            * self.unit
        )

    def to_string(self, **kwargs):
        return (
            (self.mean[0] * self.unit).to_string(**kwargs)
            if len(self) == 1
            else (self.mean * self.unit).to_string(**kwargs)
        )

    def _repr_latex_(self):
        return (
            f"{self.name} ({self.disttype}) {self.to_string(format='latex', subfmt='inline')} $\pm$ "
            + f"{(self.diff[0] * self.unit if len(self) == 1 else self.diff * self.unit).to_string(format='latex', subfmt='inline')}"
        )
