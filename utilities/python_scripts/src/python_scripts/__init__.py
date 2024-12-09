from enum import Enum, auto


class Gamma(Enum):
    ONELINK = auto()
    LOCAL = auto()
    VEC_ONELINK = auto()
    VEC_LOCAL = auto()
    PION_LOCAL = auto()


class ConfigBase:
    pass


PARALLEL_LOAD = False


def set_parallel_load(value: bool) -> None:
    globals()['PARALLEL_LOAD'] = value
