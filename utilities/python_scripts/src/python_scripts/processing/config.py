import typing as t
from dataclasses import dataclass

from python_scripts import ConfigBase, utils


@dataclass
class LoadH5Config:
    """Parameters providing index names and values to convert hdf5 datasets
    into a DataFrame.

    Properties
    ----------
        name: str
            The name to give the datasets provided in `datasets`
        datasets: dict(str, str|list(str))
            Dictionary keys will correspond to DataFrame index labels.
            Dictionary values are hdf5 file paths to access corresponding data. If given
                a list of paths, the first valid path will be used.
    """
    name: str
    datasets: t.Dict[str, t.List[str]]

    def __init__(self, **kwargs):
        self.name = kwargs['name']
        self.datasets = {}
        for k,v in kwargs['datasets'].items():
            if isinstance(v,str):
                self.datasets[k] = [v]
            else:
                assert isinstance(v, t.List)
                self.datasets[k] = v


@dataclass
class LoadArrayConfig:
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
        self.labels = utils.process_params(**self.labels)


@dataclass
class DataioConfig(ConfigBase):
    filestem: str
    array_params: t.Optional[
        t.Union[LoadArrayConfig,
                t.Dict[str, LoadArrayConfig]]
    ] = None
    replacements: t.Optional[t.Dict[str, t.Union[str, t.List[str]]]] = None
    regex: t.Optional[t.Dict[str, str]] = None
    h5_params: t.Optional[LoadH5Config] = None
    dict_labels: t.Optional[t.List[str]] = None
    actions: t.Optional[t.Dict[str, t.Any]] = None

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

        self.replacements = utils.process_params(**self.replacements)

    @classmethod
    def create(cls, **kwargs):
        """Returns an instance of DataioConfig from `params` dictionary.

                Parameters
                ----------
                kwargs
                    keys should correspond to class parameters (above).
                    `h5_params` and `array_params`, if provided,
                    should have dictionaries that can be passed to `create` static methods
                    in H5Params and ArrayParams, respectively.
                """


        obj_vars = kwargs.copy()
        h5_params = obj_vars.pop('h5_params', {})
        array_params = obj_vars.pop('array_params', {})
        if h5_params:
            obj_vars['h5_params'] = LoadH5Config(**h5_params)

            obj_vars['array_params'] = {}
            for k, v in array_params.items():
                obj_vars['array_params'][k] = LoadArrayConfig(**v)
        elif array_params:
            obj_vars['array_params'] = LoadArrayConfig(**array_params)

        return DataioConfig(**obj_vars)

def get_config_factory(config_label: str):
    configs = {
        "dataio": DataioConfig.create
    }

    if config_label in configs:
        return configs[config_label]
    else:
        raise ValueError(f"No config implementation for `{config_label}`.")

def get_dataio_config(params: t.Dict) -> DataioConfig:
    return get_config_factory('dataio')(**params)
