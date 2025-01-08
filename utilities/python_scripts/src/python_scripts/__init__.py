import logging
import typing as t
import sys
from enum import Enum, auto


class Gamma(Enum):
    ONELINK = auto()
    LOCAL = auto()
    VEC_ONELINK = auto()
    VEC_LOCAL = auto()
    PION_LOCAL = auto()

    @property
    def gamma_list(self) -> t.List[str]:
        if self in [Gamma.ONELINK, Gamma.VEC_ONELINK]:
            return ["GX_G1", "GY_G1", "GZ_G1"]
        if self == Gamma.LOCAL:
            return ["G5_G5", "GX_GX", "GY_GY", "GZ_GZ"]
        if self == Gamma.VEC_LOCAL:
            return ["GX_GX", "GY_GY", "GZ_GZ"]
        if self == Gamma.PION_LOCAL:
            return ["G5_G5"]
        else:
            raise ValueError(f"Unexpected Gamma value for mesons: {self}")

    @property
    def gamma_string(self) -> str:
        gammas = self.gamma_list
        gammas = [f"({gamma})" for gamma in gammas]
        gammas = " ".join(gammas)
        gammas = gammas.replace("_", " ")
        return gammas

    @property
    def local(self) -> bool:
        if self in [Gamma.LOCAL, Gamma.PION_LOCAL, Gamma.VEC_LOCAL]:
            return True
        else:
            return False

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
