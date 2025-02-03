import logging
import typing as t
import sys
from dataclasses import fields
from enum import Enum, auto


class Gamma(Enum):
    ONELINK = auto()
    LOCAL = auto()
    VEC_ONELINK = auto()
    VEC_LOCAL = auto()
    PION_LOCAL = auto()
    G1_G1 = auto()
    GX_GX = auto()
    GY_GY = auto()
    GZ_GZ = auto()
    G5_G5 = auto()
    GX_G1 = auto()
    GY_G1 = auto()
    GZ_G1 = auto()

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
            return [self.name]
            # raise ValueError(f"Unexpected Gamma value: {self}")

    @property
    def gamma_string(self) -> str:
        gammas = self.gamma_list
        gammas = [f"({gamma})" for gamma in gammas]
        gammas = " ".join(gammas)
        gammas = gammas.replace("_", " ")
        return gammas

    @property
    def _local_gammas(self) -> t.List:
        return [Gamma.LOCAL, Gamma.PION_LOCAL, Gamma.VEC_LOCAL,
                Gamma.G1_G1,
                Gamma.GX_GX,
                Gamma.GY_GY,
                Gamma.GZ_GZ,
                Gamma.G5_G5
                ]

    @property
    def local(self) -> bool:
        if self in self._local_gammas:
            return True
        else:
            return False


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


T = t.TypeVar('T', bound='ConfigBase')


class ConfigBase:

    @classmethod
    def create(cls: t.Type[T], **kwargs) -> T:
        """Creates a new instance of ConfigBase from a dictionary.

        Note
        ----
        Checks for dataclass fields stored in the class object.
        If kwargs match a dataclass field apart from a leading underscore,
        e.g. kwargs['mass'] and _mass, the value gets assigned to the underscored object attribute.
        """

        conflicts = [k for k in kwargs.keys() if k in kwargs and f"_{k}" in kwargs and not k.startswith('_')]

        if any(conflicts):
            raise ValueError(
                f"Conflict in parameters. Both {conflicts[0]} and _{conflicts[0]} passed to {cls} `create`.")

        class_vars = [f.name for f in fields(cls)]
        obj_vars = {}
        new_vars = {}
        for k, v in kwargs.items():
            if k in class_vars:
                obj_vars[k] = v
            elif f"_{k}" in class_vars and not k.startswith('_'):
                obj_vars[f"_{k}"] = v
            elif k in cls.__dict__:
                raise ValueError(f"Cannot overwrite existing {cls} param, `{k}`. Try relabeling `{k}`")
            else:
                new_vars[k] = v

        obj = cls(**obj_vars)

        for k, v in new_vars.items():
            setattr(obj, k, v)

        return obj

    def string_dict(self):
        """Converts all attributes without leading underscore to strings or lists of strings.
        Dictionary attributes are removed from output.
        Returns a dictionary keyed by the attribute labels
        """
        res = {}
        for k, v in self.__dict__.items():
            if k.startswith("_") or isinstance(v, t.Dict):
                continue
            elif isinstance(v, t.List):
                res[k] = list(map(str, v))
            elif isinstance(v, bool):
                res[k] = str(v).lower()
            else:
                if v is not None:
                    res[k] = str(v)

        return res

    @property
    def public_dict(self):
        """Converts object attributes to dictionary, removing attributes with leading underscore
        Returns a dictionary keyed by the attribute labels
        """
        res = {}
        for k, v in self.__dict__.items():
            if not k.startswith("_"):
                res[k] = v

        return res
