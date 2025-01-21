import importlib
from time import sleep

from python_scripts.nanny import config, TaskBase, SubmitConfig
import typing as t


def get_builder_module(job_type: str, task_type: t.Optional[str] = None):
    module_path = 'python_scripts.nanny.runio'
    if job_type != 'hadrons':
        raise NotImplementedError

    module_path += f'.{job_type}'
    if task_type:
        module_path += f'.{task_type}'

    builder = importlib.import_module(module_path)

    return builder

def build_params(job_type: str, task_type: str, submit_config: SubmitConfig,
                 task_config: TaskBase, *args, **kwargs) -> t.Tuple[t.List[t.Dict],t.Optional[t.List[str]]]:
    builder = get_builder_module(job_type, task_type)

    return builder.build_params(submit_config, task_config, *args, **kwargs)


def bad_files(job_type: str, task_type: str, submit_config: SubmitConfig, task_config: TaskBase, *args, **kwargs) -> t.List[str]:
    builder = get_builder_module(job_type,task_type)

    return builder.bad_files(submit_config, task_config, *args, **kwargs)

def get_task_factory(job_type: str, task_type: str) -> t.Callable[...,TaskBase]:
    builder = get_builder_module(job_type, task_type)

    return builder.get_task_factory()

def get_submit_factory(job_type: str) -> t.Callable[..., SubmitConfig]:
    builder = get_builder_module(job_type)

    return builder.get_submit_factory()

if __name__ == '__main__':
    from python_scripts import utils

    param = utils.load_param('params.yaml')

    jc = config.get_job_config(param, 'L')
    sc = config.get_submit_config(param,jc, series='a', cfg='100')

    fc = config.get_outfile_config(param)
    stuff = config.build_params(sc, jc,fc)

    sleep(1)

