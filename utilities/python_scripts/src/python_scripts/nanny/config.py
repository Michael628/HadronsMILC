import os.path
import typing as t
from dataclasses import dataclass, field, fields, MISSING

from python_scripts import (
    Gamma, utils, config as c
)

# SUBMISSION CLASSES
# ============Submit Configuration===========
@dataclass
class SubmitConfig(c.ConfigBase):
    ens: str
    time: int
    series: str
    cfg: str

# ============Submit Contraction Configuration===========
@c.dataclass_with_getters
class SubmitContractConfig(SubmitConfig):
    _diagram_params: t.Dict
    hardware: t.Optional[str] = None
    logging_level: t.Optional[str] = None

# ============Submit Hadrons Configuration===========
@c.dataclass_with_getters
class SubmitHadronsConfig(SubmitConfig):
    tstart: int = 0
    eigresid: float = 1e-8
    blocksize: int = 500
    multifile: bool =  False
    dt: t.Optional[int] = None
    eigs: t.Optional[int] = None
    sourceeigs: t.Optional[int] = None
    noise: t.Optional[int] = None
    tstop: t.Optional[int] = None
    alpha: t.Optional[float] = None
    beta: t.Optional[int] = None
    npoly: t.Optional[int] = None
    nstop: t.Optional[int] = None
    nk: t.Optional[int] = None
    nm: t.Optional[int] = None
    _run_id: str = ''
    _mass: t.Dict[str,float] = field(default_factory=dict)
    _overwrite_sources: bool = True

    def __post_init__(self):
        if not self.mass:
            self.mass = {}
        self._mass['zero'] = 0.0

        if self.eigs:
            if not self.sourceeigs:
                self.sourceeigs = self.eigs
            if not self.nstop:
                self.nstop = self.eigs

        if self.time and not self.tstop:
            self.tstop = self.time - 1


    @property
    def run_id(self):
        return self._run_id.format(**self.string_dict)

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, value: t.Dict[str, float]) -> t.Dict[str, float]:
        assert isinstance(value,t.Dict)
        self._mass = value
        self._mass['zero'] = 0.0

    @property
    def tsource_range(self) -> t.List[int]:
        return list(range(self.tstart,self.tstop+1,self.dt))

    @property
    def mass_out_label(self):
        res = {}
        for k, v in self.mass.items():
            res[k] = str(v).removeprefix('0.')
        return res


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
            "contract": ".{cfg}.p"
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



# TASK CLASSES
# ============Task Base===========
@dataclass
class TaskBase:
    pass

# ============Custom Task===========
@dataclass
class CustomTask(TaskBase):
    file: str

# ============LMI Task===========
@dataclass
class LMITask(TaskBase):

    # ============Epack===========
    @dataclass
    class EpackTask:
        load: bool
        multifile: bool = False
        save_eigs: bool = False
        save_evals: bool = True

    # ============Operator List===========
    @dataclass
    class OpList:
        """Configuration for a list of gamma operations.

        Attributes
        ----------
        operations: list
            Gamma operations to be performed, usually for meson fields or high mode solves.
        """

        @dataclass
        class Op:
            """Parameters for a gamma operation and associated masses.
            """
            gamma: Gamma
            mass: t.List[str]

        operations: t.List[Op]

        def __init__(self, **kwargs):
            """Creates a new instance of OpList from a dictionary.

             Note
             ----
             Valid dictionary input formats:

             kwargs = {
               'gamma': ['op1','op2','op3'],
               'mass': ['m1','m2']
             }

             or

             kwargs = {
               'op1': {
                 'mass': ['m1']
               },
               'op2': {
                 'mass': ['m2','m3']
               }
             }

            """
            if 'mass' not in kwargs:
                operations = []
                for key, val in kwargs.items():
                    mass = val['mass']
                    if isinstance(mass, str):
                        mass = [mass]
                    gamma = Gamma[key.upper()]
                    operations.append(self.Op(gamma=gamma, mass=mass))
            else:
                gammas = kwargs['gamma']
                mass = kwargs['mass']
                if isinstance(mass, str):
                    mass = [mass]
                if isinstance(gammas, str):
                    gammas = [gammas]
                operations = [
                    self.Op(gamma=Gamma[g.upper()], mass=mass)
                    for g in gammas
                ]
            self.operations = operations

        @property
        def mass(self):
            res: t.Set = set()
            for op in self.operations:
                for m in op.mass:
                    res.add(m)

            return list(res)

    @dataclass
    class HighModes(OpList):
        skip_cg: bool = False

        def __init__(self, **kwargs):
            obj_vars = kwargs.copy()

            self.skip_cg = obj_vars.pop('skip_cg',self.skip_cg)

            super().__init__(**obj_vars)

    epack: t.Optional[EpackTask] = None
    meson: t.Optional[OpList] = None
    high_modes: t.Optional[HighModes] = None

    def __init__(self, **kwargs):
        """Creates a new instance of LMITaskConfig from a dictionary.
        """
        for f in fields(self):
            field_name = f.name
            field_type = t.get_args(f.type)[0]
            if field_name in kwargs:
                setattr(self,field_name,field_type(**kwargs[field_name]))


    @property
    def mass(self):
        """Returns list of labels for masses required by task components."""
        res = []

        if self.epack and not self.epack.load:
            res.append('zero')
        if self.meson:
            res += self.meson.mass
        if self.high_modes:
            res += self.high_modes.mass

        return list(set(res))


# ============Contraction Task===========
@dataclass
class ContractTask(TaskBase):
    diagrams: t.List[str]


# JOB CLASS
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
                task_type = kwargs.get('task_type',self.task_type)
                self.tasks = get_config_factory(task_type)(**kwargs['tasks'])
            elif field_name == 'infile':
                self.infile = kwargs['io']
            else:
                setattr(self,field_name,kwargs.get(field_name,field_default))

    def __post_init__(self):
        if 'run_id' not in self.params:
            self.params['run_id'] = "LMI-RW-series-{series}-{eigs}-eigs-{noise}-noise"

    def get_infile(self, submit_config: SubmitHadronsConfig) -> str:
        ext = {
            'smear':"{series}.{cfg}.txt",
            'hadrons':"{series}.{cfg}.xml",
            'contract':"{series}.{cfg}.yaml"
        }
        return f"{self.infile}-{ext[self.job_type]}".format(**submit_config.string_dict)


# ============Convenience functions===========
def get_config_factory(config_label: str):
    configs = {
        'submit_smear': SubmitConfig.create,
        'submit_hadrons': SubmitHadronsConfig.create,
        'submit_contract': SubmitContractConfig.create,
        'job_config': JobConfig,
        'outfile': OutfileList,
        'custom': CustomTask,
        'lmi': LMITask,
        'smear': c.ConfigBase.create,
        'contract': ContractTask,
    }

    if config_label in configs:
        return configs[config_label]
    else:
        raise ValueError(f"No config implementation for `{config_label}`.")

def get_job_config(param: t.Dict, step: str) -> JobConfig:
    return get_config_factory('job_config')(**param['job_setup'][step])

def get_submit_config(param: t.Dict, job_config: JobConfig, **kwargs) -> SubmitConfig:
    submit_params = utils.deep_copy_dict(param['submit_params'])
    submit_type = 'submit_' + job_config.job_type
    additional_params = job_config.job_type + '_params'
    if additional_params in param:
        submit_params.update(param[additional_params])

    if job_config.params:
        submit_params.update(job_config.params)

    return get_config_factory(submit_type)(**submit_params, **kwargs)

def get_outfile_config(param: t.Dict):
    return get_config_factory('outfile')(**param['files'])


if __name__ == '__main__':
    a = LMITask(**{
        'epack': {'load':False},
        'meson':{'gamma':'onelink','mass':'l'}
    })

