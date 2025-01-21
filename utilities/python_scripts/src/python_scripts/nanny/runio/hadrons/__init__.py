import typing as t
from dataclasses import field

from python_scripts import config as c
from python_scripts.nanny import SubmitConfig


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
    seed: t.Optional[str] = None

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
        return self._run_id.format(**self.string_dict())

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


def get_submit_factory() -> t.Callable[..., SubmitHadronsConfig]:
    return SubmitHadronsConfig.create