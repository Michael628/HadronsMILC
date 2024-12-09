import typing as t
from dataclasses import dataclass

from python_scripts import ConfigBase, Gamma


@dataclass
class EpackTaskConfig(ConfigBase):
    load: bool
    multifile: bool = False
    save_eval: t.Optional[bool] = None


def create_epack_config(params: t.Dict) -> ConfigBase:
    return EpackTaskConfig(**params)


@dataclass
class OpConfig(ConfigBase):
    """Configuration for a list of gamma operations and associated masses.
    Usually for the sake of performing a calculation for each `gamma` at each
    `mass`
    """
    gamma: t.List[Gamma]
    mass: t.List[str]


def create_op_config(params: t.Dict) -> ConfigBase:
    m = params['mass']
    if isinstance(m, str):
        m = [m]
    gamma = [Gamma[g.upper()] for g in params['gamma']],

    return OpConfig(gamma=gamma, mass=m)


@dataclass
class OpListConfig(ConfigBase):
    """Configuration for a list of gamma operations.

    Attributes
    ----------
    gammas: list
        List of gamma operations to be run.
    """
    gammas: t.List[OpConfig]


def create_op_list_config(params: t.List) -> ConfigBase:
    """Converts a list of strings matching names in Gamma enum
    (case insensitive) to corresponding gamma operators.
    Currently assumes the use of mass=`l`, from RunConfig
    """
    return OpListConfig(gammas=[OpConfig.create(g) for g in params])
