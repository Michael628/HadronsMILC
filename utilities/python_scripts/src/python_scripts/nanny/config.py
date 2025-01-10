import os.path
import typing as t
from dataclasses import dataclass, field
from python_scripts import (
    Gamma, utils
)
from python_scripts.config import ConfigBase

# ============Run Configuration===========
@dataclass
class RunConfig(ConfigBase):
    ens: str
    series: str = ''
    cfg: str = ''
    dt: t.Optional[int] = None
    eigs: t.Optional[int] = None
    sourceeigs: t.Optional[int] = None
    noise: t.Optional[int] = None
    mass: t.Dict[str,float] = field(default_factory=dict)
    time: t.Optional[int] = None
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
    overwrite_sources: bool = True

    def __post_init__(self):
        if not self.mass:
            self.mass = {}
        self.mass['zero'] = 0.0

        if self.eigs:
            if not self.sourceeigs:
                self.sourceeigs = self.eigs
            if not self.nstop:
                self.nstop = self.eigs

        if self.time and not self.tstop:
            self.tstop = self.time - 1


    @property
    def tsource_range(self) -> t.List[int]:
        return list(range(self.tstart,self.tstop+1,self.dt))

    @property
    def mass_out_label(self):
        res = {}
        for k, v in self.mass.items():
            res[k] = str(v).removeprefix('0.')
        return res

    @property
    def string_dict(self):
        """Converts all attributes to strings or lists of strings.
        Dictionary attributes are removed from output.
        Returns a dictionary keyed by the attribute labels
        """
        res = {}
        for k, v in self.__dict__.items():
            if isinstance(v,t.Dict):
                continue
            elif isinstance(v,t.List):
                res[k] = list(map(str,v))
            elif isinstance(v,bool):
                res[k] = str(v).lower()
            else:
                res[k] = str(v)

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
    def create(cls, params: t.Dict):
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
        home = params['home']
        for k in extensions:
            if k in params:
                outfiles[k] = OutfileConfig(
                    filestem=str(os.path.join(home,params[k]['filestem'])),
                    ext=extensions[k],
                    good_size=params[k]['good_size']
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
    def create(cls, params: t.Dict):
        """Creates a new instance of OpConfig from a dictionary."""
        m = params['mass']
        if isinstance(m, str):
            m = [m]
        gamma = Gamma[params['gamma'].upper()]

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
    def create(cls, params: t.Dict):
        """Creates a new instance of OpListConfig from a dictionary.

         Note
         ----
         Valid dictionary input formats:

         params = {
           'gamma': ['op1','op2','op3'],
           'mass': ['m1','m2']
         }

         or

         params = {
           'op1': {
             'mass': ['m1']
           },
           'op2': {
             'mass': ['m2','m3']
           }
         }

        """
        if 'mass' not in params:
            operations = [
                OpConfig.create({'gamma':key,'mass':val['mass']})
                for key, val in params.items()
            ]
        else:
            gammas = params['gamma']
            if isinstance(gammas, str):
                gammas = [gammas]
            operations = [
                OpConfig.create({'gamma':gamma,'mass':params['mass']})
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
    def create(cls, params: t.Dict):
        """Creates a new instance of LMITaskConfig from a dictionary.
        """
        # Assumes Valid class attributes labeled by corresponding strings in
        # `get_config_factory` function.
        config_params = {
            key: get_config_factory(key)(val)
            for key, val in params.items()
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


# ============Job Configuration===========
@dataclass
class JobConfig(ConfigBase):
    tasks: ConfigBase
    job_type: str
    infile: str
    wall_time: str
    run: str
    run_params: t.Optional[t.Dict] = None

    @classmethod
    def create(cls, params: t.Dict):
        """Creates a new instance of JobConfig from a dictionary."""
        res = {
            'tasks': get_config_factory(params['job_type'])(params['tasks']),
        }
        res.update({k:v for k,v in params.items() if k not in res})

        # Rename `input` label from parameter file to `infile` attribute
        res['infile'] = res.pop('io')

        return JobConfig(**res)

    def get_infile(self, run_config: RunConfig) -> str:
        ext = {
            'lmi':"{series}.{cfg}.xml"
        }
        return f"{self.infile}-{ext[self.job_type]}".format(**run_config.string_dict)


# ============Convenience functions===========
def get_config_factory(config_label: str):
    configs = {
        'run_config': RunConfig.create,
        'job_config': JobConfig.create,
        "epack": EpackTaskConfig.create,
        "meson": OpListConfig.create,
        "high_modes": OpListConfig.create,
        'lmi': LMITaskConfig.create,
        'outfile': OutfileListConfig.create,
    }

    if config_label in configs:
        return configs[config_label]
    else:
        raise ValueError(f"No config implementation for `{config_label}`.")

def get_job_config(param: t.Dict, step: str) -> JobConfig:
    return get_config_factory('job_config')(param['job_setup'][step])

def get_run_config(param: t.Dict, job_config: JobConfig) -> RunConfig:
    run_params = utils.deep_copy_dict(param['run_params'])
    if job_config.run_params:
        run_params.update(job_config.run_params)

    return get_config_factory('run_config')(run_params)

def get_outfile_config(param: t.Dict):
    return get_config_factory('outfile')(param['files'])
