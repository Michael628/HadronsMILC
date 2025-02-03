import logging
import typing as t
from dataclasses import dataclass

from python_scripts import config as c
from python_scripts.nanny import SubmitConfig, TaskBase
from python_scripts.nanny.config import OutfileList


@c.dataclass_with_getters
class SubmitContractConfig(SubmitConfig):
    _diagram_params: t.Dict
    hardware: t.Optional[str] = None
    logging_level: t.Optional[str] = None


@dataclass
class ContractTask(TaskBase):
    diagrams: t.List[str]

def bad_files(submit_config: SubmitContractConfig,
              task_config: ContractTask, outfile_config_list: OutfileList) -> t.List[str]:
    logging.warning(
        f"Check completion succeeds automatically. No implementation of bad_files function in `{__file__}`.")
    return []


def get_task_factory():
    return ContractTask

def get_submit_factory() -> t.Callable[..., SubmitContractConfig]:
    return SubmitContractConfig.create