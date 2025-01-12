import os.path
import typing as t
from dataclasses import dataclass, field
from python_scripts import (
    Gamma, utils
)
from python_scripts.config import ConfigBase


# ============Submit Configuration===========
@dataclass
class SubmitConfig(ConfigBase):
    ens: str
    time: int
    series: str
    cfg: str

# ============Submit Hadrons Configuration===========
@dataclass
class SubmitHadronsConfig(SubmitConfig):
    dt: t.Optional[int] = None
    eigs: t.Optional[int] = None
    sourceeigs: t.Optional[int] = None
    noise: t.Optional[int] = None
    tstart: int = 0
    tstop: t.Optional[int] = None
    alpha: t.Optional[float] = None
    beta: t.Optional[int] = None
    npoly: t.Optional[int] = None
    nstop: t.Optional[int] = None
    nk: t.Optional[int] = None
    nm: t.Optional[int] = None
    multifile: bool =  False
    eigresid: float = 1e-8
    blocksize: int = 500
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

    @property
    def tsource_range(self) -> t.List[int]:
        return list(range(self.tstart,self.tstop+1,self.dt))

    @property
    def mass_out_label(self):
        res = {}
        for k, v in self.mass.items():
            res[k] = str(v).removeprefix('0.')
        return res


# ============Outfile Configuration===========
@dataclass
class OutfileConfig(ConfigBase):
    filestem: str
    ext: str
    good_size: int


@dataclass
class OutfileListConfig(ConfigBase):
    fat_links: t.Optional[OutfileConfig] = None
    long_links: t.Optional[OutfileConfig] = None
    gauge_links: t.Optional[OutfileConfig] = None
    eig: t.Optional[OutfileConfig] = None
    eigdir: t.Optional[OutfileConfig] = None
    eval: t.Optional[OutfileConfig] = None
    high_modes: t.Optional[OutfileConfig] = None
    meson: t.Optional[OutfileConfig] = None

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
                outfiles[k] = OutfileConfig(
                    filestem=str(os.path.join(home,kwargs[k]['filestem'])),
                    ext=extensions[k],
                    good_size=kwargs[k]['good_size']
                )

        return OutfileListConfig(**outfiles)


# ============Epack Task Configuration===========
@dataclass
class EpackTaskConfig(ConfigBase):
    load: bool
    multifile: bool = False
    save_eigs:  bool = False
    save_evals: bool = True


# ============Operator Configuration===========
@dataclass
class OpConfig(ConfigBase):
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

        return OpConfig(gamma=gamma, mass=m)


@dataclass
class OpListConfig(ConfigBase):
    """Configuration for a list of gamma operations.

    Attributes
    ----------
    operations: list
        List of gamma operations to be run.
    """
    operations: t.List[OpConfig]

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
                OpConfig.create(gamma=key,mass=val['mass'])
                for key, val in kwargs.items()
            ]
        else:
            gammas = kwargs['gamma']
            if isinstance(gammas, str):
                gammas = [gammas]
            operations = [
                OpConfig.create(gamma=gamma,mass=kwargs['mass'])
                for gamma in gammas
            ]

        return OpListConfig(operations=operations)
    pass

    @property
    def mass(self):
        res: t.Set = set()
        for op in self.operations:
            for m in op.mass:
                res.add(m)

        return list(res)



# ============LMI Task Configuration===========
@dataclass
class LMITaskConfig(ConfigBase):
    epack: t.Optional[EpackTaskConfig] = None
    meson: t.Optional[OpListConfig] = None
    high_modes: t.Optional[OpListConfig] = None


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
        return LMITaskConfig(**config_params)
    pass

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


# ============Submit Contraction Configuration===========
@dataclass
class ContractTaskConfig(ConfigBase):
    diagrams: t.List[str]

@dataclass
class SubmitContractConfig(SubmitConfig):
    diagram_params: t.Dict
    hardware: t.Optional[str] = None
    logging_level: t.Optional[str] = None

# ============Job Configuration===========
@dataclass
class JobConfig(ConfigBase):
    tasks: ConfigBase
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
        "epack": EpackTaskConfig.create,
        "meson": OpListConfig.create,
        "high_modes": OpListConfig.create,
        'lmi': LMITaskConfig.create,
        'outfile': OutfileListConfig.create,
        'smear': ConfigBase.create,
        'contract': ContractTaskConfig.create
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
