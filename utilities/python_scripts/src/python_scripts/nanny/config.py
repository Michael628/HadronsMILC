import os.path
import typing as t
from dataclasses import dataclass, field, fields, MISSING

from python_scripts import utils
from python_scripts.nanny import TaskBase, runio, SubmitConfig

# OUTFILE CLASSES
# ============Outfile Parameters===========
@dataclass
class OutfileList:

    @dataclass
    class Outfile:
        filestem: str
        ext: str
        good_size: int

    fat_links: t.Optional[Outfile] = None
    long_links: t.Optional[Outfile] = None
    gauge_links: t.Optional[Outfile] = None
    eig: t.Optional[Outfile] = None
    eigdir: t.Optional[Outfile] = None
    eval: t.Optional[Outfile] = None
    high_modes: t.Optional[Outfile] = None
    meson: t.Optional[Outfile] = None
    a2avec: t.Optional[Outfile] = None

    def __init__(self, **kwargs):
        """Creates a new instance of OutfileConfigList from a dictionary."""
        extensions = {
            "fat_links": ".{cfg}",
            "long_links": ".{cfg}",
            "gauge_links": ".{cfg}",
            "eig": ".{cfg}.bin",
            "eigdir": ".{cfg}/v{eig_index}.bin",
            "eval": ".{cfg}.h5",
            "high_modes": ".{cfg}.h5",
            "meson": ".{cfg}/{gamma}_0_0_0.h5",
            "contract": ".{cfg}.p",
            "a2avec": ".{cfg}.bin"
        }
        home = kwargs['home']
        for k in extensions:
            if k in kwargs:
                outfile = self.Outfile(
                    filestem=str(os.path.join(home,kwargs[k]['filestem'])),
                    ext=extensions[k],
                    good_size=kwargs[k]['good_size']
                )
                setattr(self,k,outfile)


# ============Job Configuration===========
@dataclass
class JobConfig:
    tasks: TaskBase
    infile: str
    wall_time: str
    run: str
    job_type: str  = 'hadrons'
    task_type: str  = 'lmi'
    params: t.Dict = field(default_factory=dict)

    def __init__(self, **kwargs):
        """Creates a new instance of JobConfig from a dictionary."""
        for f in fields(self):
            field_name = f.name

            field_default = None
            if f.default_factory is not MISSING:
                field_default = f.default_factory()
            elif f.default is not MISSING:
                field_default = f.default

            if field_name == 'tasks':
                continue
            elif field_name == 'infile':
                self.infile = kwargs['io']
            else:
                setattr(self,field_name,kwargs.get(field_name,field_default))

        task_params = kwargs.get('tasks', {})
        task_builder: t.Callable = runio.get_task_factory(self.job_type, self.task_type)
        self.tasks = task_builder(**task_params)

        if 'run_id' not in self.params:
            self.params['run_id'] = "LMI-RW-series-{series}-{eigs}-eigs-{noise}-noise"

    def get_infile(self, submit_config: SubmitConfig) -> str:
        ext = {
            'smear':"{series}.{cfg}.txt",
            'hadrons':"{series}.{cfg}.xml",
            'contract':"{series}.{cfg}.yaml"
        }
        return f"{self.infile}-{ext[self.job_type]}".format(**submit_config.string_dict())


# ============Convenience functions===========
def get_job_config(param: t.Dict, step: str) -> JobConfig:
    return JobConfig(**param['job_setup'][step])

def get_submit_config(param: t.Dict, job_config: JobConfig, **kwargs) -> SubmitConfig:
    submit_params = utils.deep_copy_dict(param['submit_params'])
    additional_params = job_config.job_type + '_params'

    if additional_params in param:
        submit_params.update(param[additional_params])
    if job_config.params:
        submit_params.update(job_config.params)

    return runio.get_submit_factory(job_config.job_type)(**submit_params, **kwargs)

def get_outfile_config(param: t.Dict):
    return OutfileList(**param['files'])

def build_params(submit_config: SubmitConfig,
                 job_config: JobConfig, *args, **kwargs) -> t.Tuple[t.List[t.Dict],t.Optional[t.List[str]]]:
    return runio.build_params(job_config.job_type,
                              job_config.task_type,
                              submit_config,
                              job_config.tasks,
                              *args, **kwargs)


def bad_files(submit_config: SubmitConfig, job_config: JobConfig, *args, **kwargs) -> t.List[str]:
    return runio.bad_files(job_config.job_type,
                           job_config.task_type,
                           submit_config,
                           job_config.tasks,
                           *args, **kwargs)


if __name__ == '__main__':
    lmi_builder = runio.get_task_factory('hadrons', 'lmi')

    a = lmi_builder(**{
        'epack': {'load':False},
        'meson':{'gamma':'onelink','mass':'l'}
    })

