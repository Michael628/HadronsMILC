import os
import pickle
import logging
import itertools
from functools import partial
import re
from string import Formatter
import numpy as np
import pandas as pd
import h5py
import python_scripts.processing.processor as processor
import typing as t


# ------ Types and config classes ------ #
class ArrayParams(t.TypedDict):
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


class H5Params(t.TypedDict):
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


class DataioConfig(t.TypedDict):
    filestem: str
    replacements: t.Dict[str, t.Union[str, t.List[str]]]
    regex: t.Dict[str, str]
    h5_params: H5Params
    array_params: t.Union[ArrayParams, t.Dict[str, ArrayParams]]
    dict_labels: t.List[str]
    actions: t.Dict[str, t.Any]


dataFrameFn = t.Callable[[np.ndarray], pd.DataFrame]
loadFn = t.Callable[[str], pd.DataFrame]
# ------ End types and config classes ------ #


# ------ IO Functions ------ #
def parse_ranges(val: t.Union[t.List, str]) -> t.List:
    """Does nothing to `val` parameters of type list, but converts
    strings to number ranges.
    """
    if isinstance(val, t.List):
        return val
    elif '..' in val:
        range_input = list(map(int, val.split("..")))
        range_input[1] += 1
        return list(range(*range_input))
    else:
        raise ValueError(
            ("`array_labels` must be lists or "
             "strings of the form `<min>..<max>`."))


def formatkeys(format_string: str) -> t.List[str]:
    """Get formatting variables found in `format_string`"""

    key_list = list(
        {
            k[1]
            for k in Formatter().parse(format_string)
            if k[1] is not None
        }
    )

    return key_list


def file_regex_gen(
        filestem: partial,
        regex: t.Dict[str, str]):
    """Formats `filestem` with replacements from `regex` and performs regex
    search on system files.

    Yields
    ------
    (replacements, filename)
        replacements: dict
            The matched values in the regex pattern corresponding to the
            formatting keys in `filestem`

        filename: str
            The file name that matched the regex search
    """
    if len(regex) == 0:
        yield {}, filestem()
    else:
        # Build regex objects to catch each replacement
        regex_repl = {k: f"(?P<{k}>{val})" for k, val in regex.items()}
        file_pattern = filestem(**regex_repl)

        # FIX ME: Assumes all regex matches occur in file name,
        # not in the directory path.
        directory, match = os.path.split(file_pattern)

        files: t.List[str] = os.listdir(directory)

        regex_pattern: re.Pattern = re.compile(match)

        for file in files:
            try:
                regex_repl = next(regex_pattern.finditer(file)).groupdict()
            except StopIteration:
                continue

            yield regex_repl, f"{directory}/{file}"


def string_replacement_gen(
        fstring: str,
        replacements: t.Dict[str, t.Union[str, t.List[str]]]):
    """Generator for keyword replacements of `fstring`

    Yields
    ------
    (repl, repl_string)
        `repl` : dict
            The replacement dictionary
            which can be passed to str.format() as kwargs
        `repl_string` : functools.partial
            The `fstring` partially formatted by `repl`. If
            no other replacements are needed then repl_string() will
            return the desired string
    """

    if len(replacements) == 0:
        yield {}, partial(fstring.format)
    else:
        keys, repls = zip(*(
            (k, map(str, r))
            if isinstance(r, t.List)
            else (k, [str(r)])
            for k, r in replacements.items()
        ))

        for r in itertools.product(*repls):
            repl: t.Dict = dict(zip(keys, r))
            string_repl: partial = partial(
                fstring.format, **repl)

            yield repl, string_repl


def ndarray_to_frame(
        array: np.ndarray,
        array_params: ArrayParams) -> pd.DataFrame:
    """Converts ndarray into a pandas DataFrame object indexed by the values in
    `array_params`. See `ArrayParams` class for details.
    """

    if len(array_params["order"]) == 1 and array_params["order"][0] == 'dt':
        array_params["labels"]['dt'] = list(range(np.prod(array.shape)))

    assert len(array_params["labels"]) == len(array_params["order"])

    indices = [array_params["labels"][k] for k in array_params["order"]]

    index: pd.MultiIndex = pd.MultiIndex.from_tuples(
        itertools.product(*indices),
        names=array_params["order"]
    )

    return pd.Series(
        array.reshape((-1,)),
        index=index,
        name='corr'
    ).to_frame()


def h5_to_frame(file: h5py.File,
                data_to_frame: t.Dict[str, dataFrameFn],
                h5_params: H5Params) -> pd.DataFrame:
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

    h5_params: H5Params
        Parameters that map keys to datasets in `file`."""
    assert all(k in data_to_frame.keys() for k in h5_params['datasets'].keys())

    df = []
    for k, v in h5_params["datasets"].items():
        frame = data_to_frame[k](file[v][:].view(np.complex128))
        if len(h5_params["datasets"]) > 1:
            frame[h5_params["name"]] = k
        df.append(frame)

    df = pd.concat(df)

    if len(h5_params["datasets"]) > 1:
        df.set_index(h5_params["name"], append=True, inplace=True)

    return df


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


def load_files(filestem: str, file_loader: loadFn,
               replacements: t.Optional[t.Dict] = None,
               regex: t.Optional[t.Dict] = None):
    """"""
    repl_keys: t.List[str] = formatkeys(filestem)

    str_repl: t.Dict = replacements if replacements else {}
    regex_repl: t.Dict = regex if regex else {}

    assert len(repl_keys) == len(str_repl) + len(regex_repl)
    assert all((
        (k in str_repl or k in regex_repl)
        for k in repl_keys
    ))

    df = []

    for str_reps, repl_filename in string_replacement_gen(
            filestem, str_repl):
        for reg_reps, regex_filename in file_regex_gen(
                repl_filename, regex_repl):
            logging.debug(f"Loading file: {regex_filename}")
            new_data: pd.DataFrame = file_loader(regex_filename)

            if len(str_reps) != 0:
                new_data[list(str_reps.keys())] = tuple(str_reps.values())
            if len(reg_reps) != 0:
                new_data[list(reg_reps.keys())] = tuple(reg_reps.values())

            df.append(new_data)

    df = pd.concat(df)

    if len(repl_keys) != 0:
        df.set_index(repl_keys, append=True, inplace=True)

    return df


def write_frame(df: pd.DataFrame, filestem: str) -> None:
    """Write DataFrame to `filestem`. If `filestem` contains format keys,
    expects columns in `df` with names matching the format keys.
    Corresponding columns will be removed from `df` and values will
    be used to format `filestem`.
    """
    if repl_keys := formatkeys(filestem):
        assert len(df) != 0
        assert all([k in df.index.names for k in repl_keys])

        for group, df_group in df.groupby(level=repl_keys):
            repl = dict(zip(repl_keys, group))

            df_group.reset_index(repl_keys, drop=True
                                 ).to_hdf(filestem.format(**repl), key='corr')

    else:
        df.to_hdf(filestem, key='corr')


def load(config: DataioConfig) -> pd.DataFrame:

    def pickle_loader(filename: str):

        dict_labels: t.Tuple = tuple(config.get("dict_labels", []))

        array_params: ArrayParams = config.get("array_params", ArrayParams())
        data_to_frame = partial(ndarray_to_frame, array_params=array_params)

        with open(filename, "rb") as f:
            data = pickle.load(f)

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
            h5_params: H5Params = config.get("h5_params")
            array_params: t.Dict[str, ArrayParams] = config.get("array_params")
            data_to_frame = {
                k: partial(ndarray_to_frame, array_params=array_params[k])
                for k in array_params.keys()
            }

            file = h5py.File(filename)

            return h5_to_frame(file, data_to_frame, h5_params)

    replacements: t.Dict = config.get("replacements", {})
    regex: t.Dict = config.get("regex", {})

    filestem: str = config.get("filestem")

    if filestem.endswith(".p"):
        file_loader = pickle_loader
    elif filestem.endswith(".h5"):
        file_loader = h5_loader
    else:
        raise ValueError("File must have extension '.p' or '.h5'")

    df = load_files(filestem, file_loader, replacements, regex)

    actions = config.get("actions", {})

    return processor.execute(df, actions=actions)

# ------ End IO functions ------ #
