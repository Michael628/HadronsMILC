import typing as t

from pydantic.dataclasses import dataclass

import python_scripts


@dataclass
class TaskBase:

    @classmethod
    def from_dict(cls, kwargs: t.Dict[str, t.Any]) -> 'TaskBase':
        return cls(**kwargs)

@dataclass
class SubmitConfig(python_scripts.ConfigBase):
    ens: str
    time: int

