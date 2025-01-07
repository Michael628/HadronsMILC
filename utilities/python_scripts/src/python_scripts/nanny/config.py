import typing as t
from dataclasses import dataclass

from python_scripts import (
    ConfigBase,
    Gamma,
)


@dataclass
class RunConfig(ConfigBase):
    ens: str
    series: str = ''
    cfg: str = ''
    dt: t.Optional[int] = None
    eigs: t.Optional[int] = None
    sourceeigs: t.Optional[int] = None
    noise: t.Optional[int] = None
    mass: t.Optional[t.Dict[str,float]] = None
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
    eigout: str = ''

    def __post_init__(self):
        if self.eigs:
            if not self.sourceeigs:
                self.sourceeigs = self.eigs
            if not self.nstop:
                self.nstop = self.eigs

        if self.time and not self.tstop:
            self.tstop = self.time - 1


    @property
    def time_range(self):
        return list(range(self.tstart,self.tstop,self.dt))

    @property
    def format_dict(self):
        """Converts all attributes to strings, lists of strings, or dictionaries of strings
        and returns a dictionary keyed by the attribute labels
        """
        res = {}
        for k, v in self.__dict__:
            if isinstance(v,t.Dict):
                res[k] = {
                    k_inner: str(v_inner)
                    for k_inner, v_inner in v.items()
                }
            elif isinstance(v,t.List):
                res[k] = list(map(str,v))
            else:
                res[k] = str(v)

        return res

@dataclass
class OutfileConfig(ConfigBase):
    filestem: str
    good_size: int


class OutfileConfigList(t.TypedDict):
    fat_links: OutfileConfig
    long_links: OutfileConfig
    gauge_links: OutfileConfig
    eig: OutfileConfig
    eigdir: OutfileConfig
    eval: OutfileConfig
    high_modes: OutfileConfig
    meson: OutfileConfig


@dataclass
class EpackTaskConfig(ConfigBase):
    load: bool
    multifile: bool = False
    save_eigs:  bool = False
    save_evals: bool = True


@dataclass
class OpConfig(ConfigBase):
    """Configuration for a list of gamma operations and associated masses.
    Usually for the sake of performing a calculation for each `gamma` at each
    `mass`
    """
    gamma: Gamma
    mass: t.List[str]

@dataclass
class OpListConfig(ConfigBase):
    """Configuration for a list of gamma operations.

    Attributes
    ----------
    gammas: list
        List of gamma operations to be run.
    """
    gammas: t.List[OpConfig]


@dataclass
class GenerateLMITaskConfig(ConfigBase):
    epack: t.Optional[EpackTaskConfig] = None
    meson: t.Optional[OpListConfig] = None
    high_modes: t.Optional[OpListConfig] = None


def get_config_factory(config_label: str):
    T = t.TypeVar('T')
    def create_config(cls: t.Type[T], params: t.Dict) -> T:
        return cls(**params)

    def create_op_config(params: t.Dict) -> OpConfig:
        m = params['mass']
        if isinstance(m, str):
            m = [m]
        gamma = Gamma[params['gamma'].upper()]

        return OpConfig(gamma=gamma, mass=m)

    def create_op_list_config(params: t.Dict) -> OpListConfig:
        """Creates a list of OpConfig objects. Supports two parameter formats.

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
            gammas = [
                create_op_config({'gamma':key,'mass':val['mass']})
                for key, val in params.items()
            ]
        else:
            gammas = [
                create_op_config({'gamma':gamma,'mass':params['mass']})
                for gamma in params['gamma']
            ]

        return OpListConfig(gammas=gammas)


    def create_generate_lmi_config(params: t.Dict) -> GenerateLMITaskConfig:
        config_params = {
            key: get_config_factory(key)(val)
            for key, val in params.items()
        }
        return GenerateLMITaskConfig(**config_params)

    def create_outfile_config(params: t.Dict) -> OutfileConfigList:
        outfiles = {}
        home = params['home']
        for k, v in params.items():
            if k != 'home':

                outfiles[k] = OutfileConfig(
                    filestem=f"{home}/{v['filestem']}",
                    good_size=v['good_size']
                )
        return OutfileConfigList(**outfiles)


    configs = {
        'run_config': lambda p: create_config(RunConfig, p),
        "epack": lambda p: create_config(EpackTaskConfig, p),
        "meson": create_op_list_config,
        "high_modes": create_op_list_config,
        'generate_lmi': create_generate_lmi_config,
        'outfile': create_outfile_config,
    }

    if config_label in configs:
        return configs[config_label]
    else:
        raise ValueError(f"No config implementation for `{config_label}`.")


def get_run_config(param: t.Dict) -> RunConfig:
    return get_config_factory('run_config')(param['lmi_param'])

def get_task_config(step: str, param: t.Dict) -> ConfigBase:
    tasks = param['job_setup'][step]['tasks']
    job_type = param['job_setup'][step]['param_file']
    return get_config_factory(job_type)(tasks)

def get_outfile_config(param: t.Dict):
    return get_config_factory('outfile')(param['files'])



def generate_outfile_formatting(task_config: ConfigBase, outfile_config: OutfileConfigList, run_config: RunConfig):
    assert isinstance(task_config,GenerateLMITaskConfig)

    if task_config.epack:
        if task_config.epack.save_eigs:
            if task_config.epack.multifile:
                yield {'eig_index': list(range(int(run_config.eigs)))}, outfile_config['eigdir']
            else:
                yield {}, outfile_config['eig']
        if task_config.epack.save_eigs:
            yield {}, outfile_config['eval']

    mass_labels = {
        key: str(val).removeprefix('0.')
        for key, val in run_config.mass.items()
    }
    res: t.Dict = {}
    for op in task_config.meson.gammas:
        if op.gamma == Gamma.ONELINK:
            gamma_label = 'G{0}_G1'
            gamma_dirs = ['X', 'Y', 'Z']
            res['gamma'] = [gamma_label.format(d) for d in gamma_dirs]
        elif op.gamma == Gamma.LOCAL:
            gamma_label = 'G{0}_G{0}'
            gamma_dirs = ['5', 'X', 'Y', 'Z']
            res['gamma'] = [gamma_label.format(d) for d in gamma_dirs]
        else:
            raise ValueError(f"Unexpected Gamma value for mesons: {op.gamma}")
        gamma_label = op.gamma.name.lower()
        res['mass'] = [mass_labels[m] for m in op.mass]
        yield res, outfile_config['meson']

    res = {'tsource': list(range(run_config.tstart, run_config.tstop, run_config.dt))}
    if task_config.epack:
        res['dset'] = ['ama','ranLL']
    else:
        res['dset'] = ['ama']

    for op in task_config.high_modes.gammas:
        res['gamma'] = op.gamma.name.lower()
        res['mass'] = [mass_labels[m] for m in op.mass]
        yield res, outfile_config['high_modes']

