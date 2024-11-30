from dataclasses import dataclass, field
import python_scripts.utils as utils
import typing as t


OutputFn = t.Callable[..., t.List[str]]


@dataclass
class ArrayParams:
    """Parameters providing index names and values to convert ndarray
    into a DataFrame.

    Properties
    ----------
        order: list
            A list of names given to each dimension of a
            multidimensional ndarray in the order of the array's shape
        labels: dict(str, Union(str, list))
            labels for each index of multidimensional array. Dictionary
            keys should match entries in `order`. Dictionary values
            should either be lists with length matching corresponding
            dimension of ndarray, or a string range in the format 'start..stop'
            where start-stop = length of array for the given dimension,
            (note: stop is inclusive).
    """
    order: t.List
    labels: t.Dict[str, t.Union[str, t.List]]

    def __post_init__(self):
        self.labels = utils.process_params(self.labels)


@dataclass
class H5Params:
    """Parameters providing index names and values to convert hdf5 datasets
    into a DataFrame.

    Properties
    ----------
        name: str
            The name to give the datasets provided in `datasets`
        datasets: dict(str, str)
            Dictionary keys will correspond to DataFrame index labels.
            Dictionary values are hdf5 file paths to access corresponding data.
    """
    name: str
    datasets: t.Dict[str, str]


@dataclass
class DataioConfig:
    filestem: str
    array_params: t.Union[ArrayParams, t.Dict[str, ArrayParams]]
    replacements: t.Optional[t.Dict[str, t.Union[str, t.List[str]]]] = None
    regex: t.Optional[t.Dict[str, str]] = None
    h5_params: t.Optional[H5Params] = None
    dict_labels: t.Optional[t.List[str]] = None
    actions: t.Optional[t.Dict[str, t.Any]] = None

    @staticmethod
    def from_dict(params: t.Dict):
        """Returns an instance of DataioConfig from `params` dictionary.

        Parameters
        ----------
        params: dict
            keys should correspond to class parameters (above).
            `h5_params` and `array_params`, if provided,
            should have dictionaries that can be passed as kwargs to H5Params
            and ArrayParams, respectively
        """

        config_params = utils.deep_copy_dict(params)

        h5_params = config_params.pop('h5_params', {})
        array_params = config_params.pop('array_params', {})
        if h5_params:
            config_params['h5_params'] = H5Params(**h5_params)

            config_params['array_params'] = {}
            for k, v in array_params:
                config_params['array_params'][k] = ArrayParams(**v)
        else:
            config_params['array_params'] = ArrayParams(**array_params)

        return DataioConfig(**config_params)

    def __post_init__(self):
        """Set defaults and process `replacements`"""

        if not self.replacements:
            self.replacements = {}
        if not self.regex:
            self.regex = {}
        if not self.dict_labels:
            self.dict_labels = []
        if not self.actions:
            self.actions = {}

        self.replacements = utils.process_params(self.replacements)


@dataclass
class Config:
    @staticmethod
    def from_dict(params: t.Dict):
        return Config(**params)


@dataclass
class EpackTaskConfig(Config):
    load: bool
    multifile: bool = False
    save_eval: t.Optional[bool] = None
    _output: t.Union[t.List[str], OutputFn] = field(init=False)

    @property
    def output(self):
        return self._output

    def __post_init__(self):
        self._output = {}

        def multifile_output(**params) -> t.List[str]:
            """Generate list of eigenvector output files from 0 to
            num eigs -1.
            """

            filestem: str = ("eigen/eig{ens}nv{eigs}er8"
                             "_grid_{series}.{cfg}/v{index}.bin")

            assert 'eigs' in params

            eigs = int(params['eigs'])
            return [filestem.format(**params, index=i) for i in range(eigs)]

        if self.multifile:
            self._output = multifile_output
        else:
            self._output['eig'] = [
                "eigen/eig{ens}nv{eigs}er8_grid_{series}.{cfg}.bin"
            ]


@dataclass
class MesonTaskConfig(Config):
    pass


@dataclass
class HighModeTaskConfig(Config):
    pass


class ConfigFactory:
    @ staticmethod
    def create_config(config_label: str, params: t.Dict):
        configs = {
            "generate_lmi": {
                "epack": EpackTaskConfig,
                "meson": MesonTaskConfig,
                "high_modes": HighModeTaskConfig
            }
        }

        if config_label in configs:
            return {
                k: v.from_dict(params)
                for k, v in configs[config_label].items()
            }
        else:
            raise ValueError(f"No config implementation for `{config_label}`.")
