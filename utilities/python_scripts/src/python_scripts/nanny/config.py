import os.path
import typing as t
from dataclasses import dataclass, field

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

    def __post_init__(self):
        assert 'zero' not in self.mass
        self.mass['zero'] = 0.0

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
        for k, v in self.__dict__:
            if isinstance(v,t.Dict):
                continue
            elif isinstance(v,t.List):
                res[k] = list(map(str,v))
            else:
                res[k] = str(v)

        return res

@dataclass
class OutfileConfig(ConfigBase):
    filestem: str
    good_size: int


@dataclass
class OutfileConfigList(ConfigBase):
    fat_links: t.Optional[OutfileConfig] = None
    long_links: t.Optional[OutfileConfig] = None
    gauge_links: t.Optional[OutfileConfig] = None
    eig: t.Optional[OutfileConfig] = None
    eigdir: t.Optional[OutfileConfig] = None
    eval: t.Optional[OutfileConfig] = None
    high_modes: t.Optional[OutfileConfig] = None
    meson: t.Optional[OutfileConfig] = None

    @property
    def eigstem(self):
        if self.eigdir:
            head, tail = os.path.split(self.eigdir.filestem)
            if tail.startswith('v') and tail.endswith('.bin'):
                return head
            else:
                return self.eigdir
        elif self.eig:
            head, tail = os.path.split(self.eig.filestem)
            tail = tail.split('.')[0]
            return os.path.join(head,tail)
        else:
            return None

    @property
    def evalstem(self):
        if self.eval:
            head, tail = os.path.split(self.eval.filestem)
            tail = tail.split('.')[0]
            return os.path.join(head,tail)
        else:
            return None

    @property
    def mesonstem(self):
        if self.meson:
            head, tail = os.path.split(self.meson.filestem)
            head, tail = os.path.split(head)
            tail = tail.split('.')[0]
            return os.path.join(head, tail)
        else:
            return None

    @property
    def highstem(self):
        if self.high_modes:
            head, tail = os.path.split(self.high_modes.filestem)
            tail = tail.split('.')[0]
            return os.path.join(head, tail)
        else:
            return None


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
    operations: list
        List of gamma operations to be run.
    """
    operations: t.List[OpConfig]


    @property
    def mass(self):
        res: t.Set = set()
        for op in self.operations:
            for m in op.mass:
                res.add(m)

        return list(res)


@dataclass
class GenerateLMITaskConfig(ConfigBase):
    epack: t.Optional[EpackTaskConfig] = None
    meson: t.Optional[OpListConfig] = None
    high_modes: t.Optional[OpListConfig] = None

    @property
    def mass(self):
        res = []

        if self.epack and not self.epack.load:
            res.append('zero')
        if self.meson:
            res += self.meson.mass
        if self.high_modes:
            res += self.high_modes.mass

        return res


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
            operations = [
                create_op_config({'gamma':key,'mass':val['mass']})
                for key, val in params.items()
            ]
        else:
            operations = [
                create_op_config({'gamma':gamma,'mass':params['mass']})
                for gamma in params['gamma']
            ]

        return OpListConfig(operations=operations)


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
                    filestem=str(os.path.join(home,v['filestem'])),
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
