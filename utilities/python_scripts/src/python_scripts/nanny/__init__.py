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

