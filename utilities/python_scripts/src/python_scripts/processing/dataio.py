import sys
import os
import logging
import itertools
from functools import partial
import numpy as np
import pandas as pd
import h5py
import typing as t

from python_scripts import utils
from python_scripts.processing import (
    processor,
    config as proc_conf
)


dataFrameFn = t.Callable[[np.ndarray], pd.DataFrame]
loadFn = t.Callable[[str], pd.DataFrame]


# ------ Data structure Functions ------ #
def ndarray_to_frame(
        array: np.ndarray,
        array_params: proc_conf.LoadArrayConfig) -> pd.DataFrame:
    """Converts ndarray into a pandas DataFrame object indexed by the values in
    `array_params`. See `ArrayParams` class for details.
    """

    if len(array_params.order) == 1 and array_params.order[0] == 'dt':
        array_params.labels['dt'] = list(range(np.prod(array.shape)))

    assert len(array_params.labels) == len(array_params.order)

    indices = [array_params.labels[k] for k in array_params.order]

    index: pd.MultiIndex = pd.MultiIndex.from_tuples(
        itertools.product(*indices),
        names=array_params.order
    )

    return pd.Series(
        array.reshape((-1,)),
        index=index,
        name='corr'
    ).to_frame()


def h5_to_frame(file: h5py.File,
                data_to_frame: t.Dict[str, dataFrameFn],
                h5_params: proc_conf.LoadH5Config) -> pd.DataFrame:
    """Converts hdf5 format `file` to pandas DataFrame based on dataset info
    provided by `h5_params`. See H5Params class for details.

    Parameters
    ----------
    file: h5py.File
        hdf5 file, the contents of which will be converted into a pandas
        DataFrame.

    data_to_frame: Dict[str, dataFrameFn]
        A dictionary of functions that convert ndarrays to DataFrames (usually,
        these are partially evaluated instances of the ndarray_to_frame
        function). The dictionary keys should match the keys of
        `h5_params['datasetes'].

    h5_params: LoadH5Config
        Parameters that map keys to datasets in `file`."""
    assert all(k in data_to_frame.keys() for k in h5_params.datasets.keys())

    df = []
    for k, v in h5_params.datasets.items():
        frame = data_to_frame[k](file[v][:].view(np.complex128))
        if len(h5_params.datasets) > 1:
            frame[h5_params.name] = k
        df.append(frame)

    df = pd.concat(df)

    if len(h5_params.datasets) > 1:
        df.set_index(h5_params.name, append=True, inplace=True)

    return df


def frame_to_dict(df: pd.DataFrame, dict_depth: int) \
        -> t.Union[t.Dict, np.ndarray]:
    num_indices = len(df.index.names)
    assert dict_depth >= 0
    assert dict_depth <= num_indices

    shape = [
        len(df.index.get_level_values(i).drop_duplicates())
        for i in range(dict_depth, num_indices)
    ]
    shape = tuple([-1] + shape) if dict_depth != 0 else tuple(shape)

    keys = [
        df.sort_index().index.get_level_values(i).drop_duplicates().to_list()
        for i in range(dict_depth)
    ]

    def join_str_fn(x):
        return ".".join(map(str, x))

    keys = list(map(join_str_fn, list(itertools.product(*keys))))

    array = df.sort_index()['corr'].to_numpy().reshape(shape)

    if dict_depth == 0:
        return array
    else:
        return {k: array[i] for k, i in zip(keys, range(len(array)))}


def dict_to_frame(
        data: t.Union[t.Dict, np.ndarray],
        data_to_frame: dataFrameFn,
        dict_labels: t.Tuple[str] = ()) -> pd.DataFrame:
    """Converts nested dictionary `data` to pandas DataFrame.

    Parameters
    ----------
    data: Union[dict, ndarray]
        Possibly nested dictionary in which leaves of dictionary tree
        are ndarrays

    data_to_frame: dataFrameFn
        A function that converts ndarrays to DataFrames (usually,
        this is a partially evaluated instance of the ndarray_to_frame
        function)

    dict_labels: tuple(str)
        Labels for each nested dictionary layer. Labels will become names of
        index levels in DataFrame and keys will become index values
        for each entry.
    """

    def entry_gen(nested: t.Dict, _index: t.Tuple = ()) \
            -> t.Generator[t.Tuple[t.Tuple, np.ndarray], None, None]:
        """Recursive Depth first search of nested dictionaries building
        list of indices from dictionary keys.


        Parameters
        ----------
            nested: dict
                The current sub-dictionary from traversing path `_index`
            _index: tuple(str)
                The sequence of keys traversed thus far in the original
                dictionary

        Yields
        ------
        (path, data)
            path: tuple(str)
                The sequence of keys traversed to get to `data` in
                the nested dictionary.
            data: ndarray
                The data that was found in `nested` by traversing indices
                in `path`.
        """

        if isinstance(next(iter(nested.values())),
                      np.ndarray):
            assert all((
                isinstance(n, np.ndarray)
                for n in nested.values()
            ))

            for key, val in nested.items():
                yield (_index + (key,), val)
        else:
            for key in nested.keys():
                yield from entry_gen(nested[key], _index + (key,))

    if isinstance(data, np.ndarray):
        return data_to_frame(data)
    else:
        assert isinstance(data, t.Dict)

        indices, concat_data = zip(*(
            (index, array) for index, array in entry_gen(data)
        ))
        concat_data = [data_to_frame(x) for x in concat_data]

        for index, frame in zip(indices, concat_data):
            frame[list(dict_labels)] = list(index)

        df = pd.concat(concat_data)

        df.set_index(list(dict_labels), append=True, inplace=True)

    return df
# ------ End data structure functions ------ #


# ------ Input functions ------ #
def load_files(filestem: str, file_loader: loadFn,
               replacements: t.Optional[t.Dict] = None,
               regex: t.Optional[t.Dict] = None):

    def proc(filename: str, repl: t.Dict) -> pd.DataFrame:
        logging.debug(f"Loading file: {filename}")
        new_data: pd.DataFrame = file_loader(filename)

        if len(repl) != 0:
            new_data[list(repl.keys())] = tuple(repl.values())

        return new_data

    return utils.process_files(filestem, proc, replacements, regex)


def load(config: proc_conf.DataioConfig) -> pd.DataFrame:

    def pickle_loader(filename: str):

        dict_labels: t.Tuple = tuple(config.dict_labels)

        array_params: proc_conf.LoadArrayConfig = config.array_params

        data_to_frame = partial(ndarray_to_frame, array_params=array_params)

        data = np.load(filename, allow_pickle=True)
        if isinstance(data, np.ndarray) and len(data.shape) == 0:
            data = data.item()

        if isinstance(data, pd.DataFrame):
            return data
        elif isinstance(data, t.Dict):
            return dict_to_frame(data,
                                 data_to_frame=data_to_frame,
                                 dict_labels=dict_labels)
        else:
            raise ValueError(
                (f"Contents of {filestem} is of type {type(data)}."
                 "Expecting dictionary or pandas DataFrame."))

    def h5_loader(filename: str):
        try:
            return pd.read_hdf(filename)
        except ValueError:
            assert config.h5_params is not None

            h5_params: proc_conf.LoadH5Config = config.h5_params

            array_params: t.Dict[str, proc_conf.LoadArrayConfig]
            array_params = config.array_params

            data_to_frame = {
                k: partial(ndarray_to_frame, array_params=array_params[k])
                for k in array_params.keys()
            }

            file = h5py.File(filename)

            return h5_to_frame(file, data_to_frame, h5_params)

    replacements: t.Dict = config.replacements
    regex: t.Dict = config.regex
    filestem: str = config.filestem

    if filestem.endswith(".p") or filestem.endswith(".npy"):
        file_loader = pickle_loader
    elif filestem.endswith(".h5"):
        file_loader = h5_loader
    else:
        raise ValueError("File must have extension '.p' or '.h5'")

    df: t.List[pd.DataFrame] = load_files(
        filestem, file_loader, replacements, regex
    )

    actions: t.Dict = config.actions
    _ = [
        processor.execute(elem, actions=actions)
        for elem in df
    ]

    df = pd.concat(df)

    return df
# ------ End Input functions ------ #


# ------ Input functions ------ #
def write_data(df: pd.DataFrame, filestem: str,
               write_fn: t.Callable[[pd.DataFrame, str], None]) -> None:
    """Write DataFrame to `filestem`. If `filestem` contains format keys,
    expects columns in `df` with names matching the format keys.
    Corresponding columns will be removed from `df` and values will
    be used to format `filestem`.

    Parameters
    ----------
    write_fn: Callable[df, filename]
        The function used to write `df` into the desired format
    """
    repl_keys = utils.formatkeys(filestem)
    if repl_keys:
        assert len(df) != 0
        assert all([k in df.index.names for k in repl_keys])

        for group, df_group in df.groupby(level=repl_keys):
            repl_vals = (group,) if isinstance(group, str) else group
            repl = dict(zip(repl_keys, repl_vals))

            filename = filestem.format(**repl)
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            write_fn(
                df_group.reset_index(repl_keys, drop=True),
                filename
            )

    else:
        write_fn(df, filestem)


def write_dict(df: pd.DataFrame, filestem: str, dict_depth: int) -> None:
    """Convert DataFrame `df` to dictionary and write to file
    based on `filestem` (see `write_data` function for details).

    Parameters
    ----------
    dict_depth: int
        If equal to 0, returns multidim ndarray with dimensions in order of df
        index levels. For `dict_depth` > 0, index levels from 0 to `dict_depth`
        are converted to dictionary keys. The remaining levels are made into
        multidim arrays.
    """
    def writeconvert(data, fname):
        np.save(fname, frame_to_dict(data, dict_depth))

    write_data(df, filestem, write_fn=writeconvert)


def write_frame(df: pd.DataFrame, filestem: str) -> None:
    """Write DataFrame `df` to hdf5 file based on `filestem`
    (see `write_data` function for details).
    """
    write_data(df, filestem,
               write_fn=lambda data, fname: data.to_hdf(fname, key='corr'))


def main(**kwargs):

    globals()['PARALLEL_LOAD'] = False
    logging_level: str
    if kwargs:
        logging_level = kwargs.pop('logging_level', 'INFO')
        config = proc_conf.get_config('load_files')(kwargs)
    else:
        try:
            params = utils.load_param('params.yaml')['load_files']
        except KeyError:
            raise ValueError("Expecting `load_files` key in params.yaml file.")

        logging_level = params.pop('logging_level', 'INFO')
        config = proc_conf.get_config('load_files')(params)

    logging.basicConfig(
        format="%(asctime)s - %(levelname)-5s - %(message)s",
        style="%",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging_level,
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
    return load(config)


if __name__ == '__main__':
    print("Assuming python interpreter is being run in interactive mode.")
    print(("Result will be stored in `result` variable"
           " upon load file completion."))
    result = main()
    logging.info("Result of file load is now stored in `result` variable.")

# ------ End Output functions ------ #
