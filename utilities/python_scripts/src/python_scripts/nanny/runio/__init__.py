import importlib
from time import sleep

from python_scripts.nanny import config
import typing as t


def get_builder_module(submit_config: config.SubmitConfig, job_config: config.JobConfig):
    module_path = 'python_scripts.nanny.runio'
    if isinstance(submit_config,config.SubmitHadronsConfig):
        module_path += '.hadrons'
        if isinstance(job_config.tasks, config.CustomTask):
            module_path += f'.{job_config.tasks.file}'
        elif isinstance(job_config.tasks, config.LMITask):
                module_path += f'.lmi'
    else:
        raise NotImplementedError

    builder = importlib.import_module(module_path)

    return builder

def build_params(submit_config: config.SubmitConfig,
                 job_config: config.JobConfig, *args, **kwargs) -> t.Tuple[t.List[t.Dict],t.Optional[t.List[str]]]:
    builder = get_builder_module(submit_config,job_config)

    return builder.build_params(submit_config, job_config.tasks, *args, **kwargs)


def bad_files(submit_config: config.SubmitConfig, job_config: config.JobConfig, *args, **kwargs) -> t.List[str]:
    builder = get_builder_module(submit_config,job_config)

    return builder.bad_files(submit_config, job_config.tasks, *args, **kwargs)

if __name__ == '__main__':
    from python_scripts import utils

    param = utils.load_param('params.yaml')

    jc = config.get_job_config(param, 'L')
    sc = config.get_submit_config(param,jc, series='a', cfg='100')
    fc = config.get_outfile_config(param)
    stuff = build_params(sc, jc,fc)

    sleep(1)