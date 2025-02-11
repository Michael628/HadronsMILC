from dataclasses import dataclass
import typing as t
import python_scripts


@dataclass
class TaskBase:
    pass


@dataclass
class SubmitConfig(python_scripts.ConfigBase):
    ens: str
    time: int
    series: str
    cfg: str
    noise: t.Optional[int] = None
    dt: t.Optional[int] = None
    eigs: t.Optional[int] = None

