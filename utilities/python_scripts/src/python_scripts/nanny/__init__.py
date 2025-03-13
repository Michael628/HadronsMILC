import typing as t

from pydantic.dataclasses import dataclass

import python_scripts


@dataclass
class TaskBase:
    pass


@dataclass
class SubmitConfig(python_scripts.ConfigBase):
    ens: str
    time: int

