import logging
import sys
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


def setup():
    logging.basicConfig(
        format="%(asctime)s - %(levelname)-5s - %(message)s",
        style="%",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
