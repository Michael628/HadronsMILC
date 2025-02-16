import typing as t
from enum import Enum, auto

import python_scripts
from python_scripts import config as c
from python_scripts.nanny.config import OutfileList

try:
    from mpi4py import MPI
    COMM = MPI.COMM_WORLD
except ImportError:
    pass


class Diagrams(Enum):
    photex = auto()
    selfen = auto()

@c.dataclass_with_getters
class DiagramConfig(python_scripts.ConfigBase):
    diagram_label: str
    contraction_type: str
    gammas: t.List[str]
    mass: str
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
    _has_high: bool = False
    _has_low: bool = False
    _npoint: int = -1
    evalfile: t.Optional[str] = None
    outfile: t.Optional[str] = None
    mesonfiles: t.Optional[str] = None


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

    @classmethod
    def create(cls, outfile_config: t.Optional[OutfileList] = None, **kwargs):
        obj_vars = kwargs.copy()
        mesonfiles = obj_vars.pop('mesonfiles')
        outfile = obj_vars.pop('outfile')
        evalfile = obj_vars.pop('evalfile',None)
        run_vars = obj_vars.pop('run_vars', {})

        def get_filename(s: str):
            if outfile_config and hasattr(outfile_config, s):
                return getattr(outfile_config, s).filename
            else:
                return s

        obj = super().create(**obj_vars)
        if isinstance(mesonfiles,str):
            mesonfiles = [mesonfiles]

        assert isinstance(mesonfiles,t.List)

        if len(mesonfiles) == 1:
            mesonfiles = mesonfiles*obj.npoint
        assert obj.npoint == len(mesonfiles)

        obj.mesonfiles = [get_filename(m) for m in mesonfiles]

        if evalfile:
            assert isinstance(evalfile,str)
            obj.evalfile = get_filename(evalfile)

        assert isinstance(outfile,str)
        obj.outfile = get_filename(outfile)

        return obj

    @property
    def meson_params(self):
        return self._meson_params

    def format_evalfile(self, **kwargs) -> None:
        self.meson_params['evalfile'] = self.evalfile.format(**kwargs)

@c.dataclass_with_getters
class RunContractConfig(python_scripts.ConfigBase):
    time: int
    ens: str
    series: str
    cfg: str
    _diagrams: t.Optional[t.List[DiagramConfig]] = None
    _overwrite_correlators: bool = True
    _hardware: str = 'cpu'
    _logging_level: str = 'INFO'
    _comm_size: int = -1
    _rank: int = -1

    @classmethod
    def create(cls, **kwargs):
        obj_vars = kwargs.copy()
        diagrams = obj_vars.pop('diagrams')
        obj = super().create(**obj_vars)
        obj.diagrams = [
            DiagramConfig.create(**v,run_vars=obj.string_dict())
            for k, v in diagrams.items()
        ]
        obj.rank = COMM.Get_rank()
        obj.comm_size = COMM.Get_size()

        return obj


def get_contract_config(params: t.Dict) -> RunContractConfig:
    return RunContractConfig.create(**params)
