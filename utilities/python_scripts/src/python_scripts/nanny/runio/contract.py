import typing as t
from dataclasses import dataclass

from python_scripts import config as c
from python_scripts.nanny import SubmitConfig, TaskBase


@c.dataclass_with_getters
class SubmitContractConfig(SubmitConfig):
    _diagram_params: t.Dict
    hardware: t.Optional[str] = None
    logging_level: t.Optional[str] = None


@dataclass
class ContractTask(TaskBase):
    diagrams: t.List[str]
