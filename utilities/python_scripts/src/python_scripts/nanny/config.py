import os.path
import typing as t
from abc import abstractmethod, ABC
from dataclasses import dataclass, field

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
class Outfile:
    filestem: str
    ext: str
    good_size: int


@dataclass
class OutfileList:
    fat_links: t.Optional[Outfile] = None
    long_links: t.Optional[Outfile] = None
    gauge_links: t.Optional[Outfile] = None
    eig: t.Optional[Outfile] = None
    eigdir: t.Optional[Outfile] = None
    eval: t.Optional[Outfile] = None
    high_modes: t.Optional[Outfile] = None
    meson: t.Optional[Outfile] = None

    @classmethod
    def create(cls, **kwargs):
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
        outfiles = {}
        home = kwargs['home']
        for k in extensions:
            if k in kwargs:
                outfiles[k] = Outfile(
                    filestem=str(os.path.join(home,kwargs[k]['filestem'])),
                    ext=extensions[k],
                    good_size=kwargs[k]['good_size']
                )

        return OutfileList(**outfiles)


# TASK CLASSES
# ============Task Base===========
class TaskBase:
    @classmethod
    def create(cls, **kwargs):
        return cls(**kwargs)


# ============Epack Task===========
@dataclass
class EpackTask(TaskBase):
    load: bool
    multifile: bool = False
    save_eigs:  bool = False
    save_evals: bool = True


# ============Operator Task===========
@dataclass
class OpTask(TaskBase):
    """Configuration for a list of gamma operations and associated masses.
    Usually for the sake of performing a calculation for each `gamma` at each
    `mass`
    """
    gamma: Gamma
    mass: t.List[str]


    @classmethod
    def create(cls, **kwargs):
        """Creates a new instance of OpConfig from a dictionary."""
        m = kwargs['mass']
        if isinstance(m, str):
            m = [m]
        gamma = Gamma[kwargs['gamma'].upper()]

        return OpTask(gamma=gamma, mass=m)


# ============Operator List Task===========
@dataclass
class OpListTask(TaskBase):
    """Configuration for a list of gamma operations.

    Attributes
    ----------
    operations: list
        List of gamma operations to be run.
    """
    operations: t.List[OpTask]

    @classmethod
    def create(cls, **kwargs):
        """Creates a new instance of OpListConfig from a dictionary.

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
            operations = [
                OpTask.create(gamma=key, mass=val['mass'])
                for key, val in kwargs.items()
            ]
        else:
            gammas = kwargs['gamma']
            if isinstance(gammas, str):
                gammas = [gammas]
            operations = [
                OpTask.create(gamma=gamma, mass=kwargs['mass'])
                for gamma in gammas
            ]

        return OpListTask(operations=operations)

    @property
    def mass(self):
        res: t.Set = set()
        for op in self.operations:
            for m in op.mass:
                res.add(m)

        return list(res)


# ============LMI Task===========
@dataclass
class LMITask(TaskBase):
    epack: t.Optional[EpackTask] = None
    meson: t.Optional[OpListTask] = None
    high_modes: t.Optional[OpListTask] = None

    @classmethod
    def create(cls, **kwargs):
        """Creates a new instance of LMITaskConfig from a dictionary.
        """
        # Assumes Valid class attributes labeled by corresponding strings in
        # `get_config_factory` function.
        config_params = {
            key: get_config_factory(key)(**val)
            for key, val in kwargs.items()
        }
        return LMITask(**config_params)

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
    job_type: str
    infile: str
    wall_time: str
    run: str
    params: t.Optional[t.Dict] = None

    @classmethod
    def create(cls, **kwargs):
        """Creates a new instance of JobConfig from a dictionary."""
        tasks = get_config_factory(kwargs['job_type'])(**kwargs['tasks'])
        res = {k:v for k,v in kwargs.items() if k not in ['tasks', 'io']}

        return JobConfig(tasks=tasks, infile=kwargs['io'], **res)

    def get_infile(self, submit_config: SubmitHadronsConfig) -> str:
        ext = {
            'smear':"{series}.{cfg}.txt",
            'lmi':"{series}.{cfg}.xml",
            'contract':"{series}.{cfg}.yaml"
        }
        return f"{self.infile}-{ext[self.job_type]}".format(**submit_config.string_dict)


# ============Convenience functions===========
def get_config_factory(config_label: str):
    configs = {
        'submit_smear': SubmitConfig.create,
        'submit_lmi': SubmitHadronsConfig.create,
        'submit_contract': SubmitContractConfig.create,
        'job_config': JobConfig.create,
        "epack": EpackTask.create,
        "meson": OpListTask.create,
        "high_modes": OpListTask.create,
        'lmi': LMITask.create,
        'outfile': OutfileList.create,
        'smear': c.ConfigBase.create,
        'contract': ContractTask.create
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
    if job_config.job_type == 'smear':
        pass
    elif job_config.job_type == 'contract':
        submit_params.update(param['contract_params'])
    else:
        submit_params.update(param['hadrons_params'])

    if job_config.params:
        submit_params.update(job_config.params)

    return get_config_factory(submit_type)(**submit_params, **kwargs)

def get_outfile_config(param: t.Dict):
    return get_config_factory('outfile')(**param['files'])
