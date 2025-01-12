import functools
import typing as t
from dataclasses import dataclass
from enum import Enum, auto
from sympy.categories import Diagram

from python_scripts.config import ConfigBase
from mpi4py import MPI

COMM = MPI.COMM_WORLD

class Diagrams(Enum):
    photex = auto()
    selfen = auto()

@dataclass
class DiagramConfig(ConfigBase):
    diagram_label: str
    contraction_type: str
    gammas: t.List[str]
    mass: str
    _outfile: functools.partial
    _mesonfile: functools.partial
    symmetric: bool = False
    newmass: t.Optional[str] = None
    high_count: t.Optional[int] = None
    high_label: t.Optional[str] = None
    low_max: t.Optional[int] = None
    low_label: t.Optional[str] = None
    mesonKey: t.Optional[str] = None
    emseedstring: t.Optional[str] = None
    perms: t.Optional[t.List[str]] = None
    n_em: t.Optional[int] = None
    _evalfile: t.Optional[str] = None
    _has_high: bool = False
    _has_low: bool = False
    _npoint: int = -1


    def __post_init__(self):

        if self.high_label and self.high_count:
            self._has_high = True

        if self.low_label:
            self._has_low = True

        npoint = {
            'conn_2pt': 2,
            'sib_conn_3pt': 3,
            'qed_conn_4pt': 4
        }
        self._npoint = npoint[self.contraction_type]

        self._meson_params = {
            "wmax_index": slice(self.low_max),
            "vmax_index": slice(self.low_max),
            "milc_mass": True
        }

        if self.newmass:
            self._meson_params['shift_mass'] = True
            self._meson_params['oldmass'] = float(f"0.{self.mass}")
            self._meson_params['newmass'] = float(f"0.{self.newmass}")
            self._meson_params['evalfile'] = self.evalfile

    @classmethod
    def create(cls, **kwargs):
        class_vars = cls.__dataclass_fields__.keys()
        public_vars = [x for x in class_vars if not x.startswith('_')]
        public_vars = {
            k:v
            for k,v in kwargs.items()
            if k in public_vars
        }
        run_vars = kwargs.get('run_vars',{})
        mesonfile = functools.partial(kwargs['mesonfile'].format,**run_vars, **public_vars)
        outfile = functools.partial(kwargs['outfile'].format, **run_vars, **public_vars)
        if 'evalfile' in kwargs:
            public_vars['_evalfile'] = kwargs['evalfile'].format(**run_vars, **public_vars)

        return DiagramConfig(_mesonfile=mesonfile,
                             _outfile=outfile,
                             **public_vars)

    @property
    def npoint(self):
        return self._npoint

    @property
    def has_high(self):
        return self._has_high

    @property
    def has_low(self):
        return self._has_low

    @property
    def meson_params(self):
        return self._meson_params

    def mesonfile(self, **kwargs) -> str:
        return self._mesonfile(**kwargs)


    def outfile(self, **kwargs) -> str:
        return self._outfile(**kwargs)


@dataclass
class RunContractConfig(ConfigBase):
    time: int
    ens: str
    series: str
    cfg: str
    _diagrams: t.List[DiagramConfig]
    _overwrite_correlators: bool = True
    _hardware: str = 'cpu'
    _logging_level: str = 'INFO'

    def __post_init__(self):
        self._rank = COMM.Get_rank()
        self._comm_size = COMM.Get_size()

    @property
    def rank(self):
        return self._rank

    @property
    def logging_level(self):
        return self._logging_level

    @property
    def hardware(self):
        return self._hardware

    @property
    def comm_size(self):
        return self._comm_size

    @property
    def diagrams(self):
        return self._diagrams

    @property
    def overwrite(self):
        return self._overwrite_correlators

    @classmethod
    def create(cls, **kwargs):
        class_vars = cls.__dataclass_fields__.keys()
        public_vars = [x for x in class_vars if not x.startswith('_')]
        private_vars = [x for x in class_vars if x.startswith('_')]

        public_vars = {
            k: v
            for k, v in kwargs.items()
            if k in public_vars
        }
        private_vars = {
            f'_{k}': v
            for k, v in kwargs.items()
            if f"_{k}" in private_vars
        }
        private_vars['_diagrams'] = [
            DiagramConfig.create(**v,run_vars=public_vars)
            for k, v in kwargs['diagrams'].items()
        ]
        return RunContractConfig(**public_vars, **private_vars)


def get_contract_config(params: t.Dict) -> RunContractConfig:
    return RunContractConfig.create(**params)
