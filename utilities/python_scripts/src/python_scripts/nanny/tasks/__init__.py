import importlib
from time import sleep

from python_scripts.nanny import config, TaskBase, SubmitConfig
import typing as t


def task_builder_module(job_type: str, task_type: t.Optional[str] = None):
    module_path = 'python_scripts.nanny.tasks'
    module_path += f'.{job_type}'
    if task_type:
        module_path += f'.{task_type}'

    builder = importlib.import_module(module_path)

    return builder

def input_params(job_type: str, task_type: str, *args, **kwargs) -> t.Tuple[t.List[t.Dict],t.Optional[t.List[str]]]:
    return task_builder_module(job_type, task_type).input_params(*args, **kwargs)

def processing_params(job_type: str, task_type: str, *args, **kwargs) -> t.Dict:
    return task_builder_module(job_type, task_type).processing_params(*args, **kwargs)

def bad_files(job_type: str, task_type: str, *args, **kwargs) -> t.List[str]:
    return task_builder_module(job_type,task_type).bad_files(*args, **kwargs)

def get_task_factory(job_type: str, task_type: str) -> t.Callable[...,TaskBase]:
    return task_builder_module(job_type, task_type).get_task_factory()

def get_submit_factory(job_type: str) -> t.Callable[..., SubmitConfig]:
    return task_builder_module(job_type).get_submit_factory()

if __name__ == '__main__':
    from python_scripts import utils

    param = utils.load_param('params.yaml')

    jc = config.get_job_config(param, 'SIB')
    sc = config.get_submit_config(param,jc, series='a', cfg='100')

    fc = config.get_outfile_config(param)
    stuff = config.input_params(jc, fc)