import os
import pickle
import logging
import itertools
import functools
import re
from string import Formatter
import numpy as np
import pandas as pd
import h5py
import python_scripts.processing.processor as processor
import typing as t
from python_scripts.utils import deep_copy_dict

DF_INTEGER_KEYS = ['cfg', 'time', 'dt']


class ArrayParams(t.TypedDict):
    order: t.List
    labels: t.Dict[str,t.Union[str,t.List]]


class HdfParams(t.TypedDict):
    name: str
    datasets: t.Dict[str,t.List[str]]


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


def formatkeys(
        format_string: str,
        keysort: t.Optional[t.Callable] = None) -> t.List[str]:

    key_list = list(
        {
            k[1]
            for k in Formatter().parse(format_string)
            if k[1] is not None
        }
    )

    if keysort is not None:
        key_list.sort(key=keysort)

    return key_list


def file_regex_gen(
        filestem: functools.partial,
        regex: t.Dict[str, str]):

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
        `repl_string` : str
            The `repl_string` string formatted by `replacements`
    """

    if len(replacements) == 0:
        yield {}, functools.partial(fstring.format)
    else:
        keys, repls = zip(*(
            (k, map(str, r))
            if isinstance(r, t.List)
            else (k, [str(r)])
            for k, r in replacements.items()
        ))

        for r in itertools.product(*repls):
            repl: t.Dict = dict(zip(keys, r))
            string_repl: functools.partial = functools.partial(
                fstring.format, **repl)

            yield repl, string_repl


def ndarray_to_frame(
        array: np.ndarray,
        array_params: ArrayParams) -> pd.DataFrame:
    """Converts ndarray into a pandas.Series object indexed by the values of
    `array_labels`, with indices named by the keys of `array_labels.
    `label_order` provides the ordering of the labels in the multidimensional
    ndarray.
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


def h5_to_df(file: h5py.File,
             h5_params: HdfParams,
             array_params: t.Dict[str,ArrayParams]) -> pd.DataFrame:
    df = []
    for k,v in h5_params["datasets"].items():
        frame = ndarray_to_frame(file[v][:].view(np.complex128), array_params[k])
        if len(h5_params["datasets"]) > 1:
            frame[h5_params["name"]] = k
        df.append(frame)

    df = pd.concat(df)

    if len(h5_params["datasets"]) > 1:
        df.set_index(h5_params["name"], append=True, inplace=True)

    return df

def dict_to_df(
        data: t.Union[t.Dict, np.ndarray],
        data_order: t.Tuple[str] = (),
        **kwargs) -> pd.DataFrame:

    def entry_gen(nested: t.Dict, _index: t.Tuple = (),
        preprocess: t.Optional[t.Callable] = None) \
            -> t.Generator[t.Tuple[t.Tuple, np.ndarray], None, None]:
        """Depth first search of nested dictionaries building
        list of indices from dictionary keys.
        """

        if isinstance(next(iter(nested.values())),
                      np.ndarray):
            assert all((
                isinstance(n, np.ndarray)
                for n in nested.values()
            ))

            for key, val in nested.items():
                yield (
                    _index+(key,),
                    preprocess(val) if preprocess else val
                )
        else:
            for key in nested.keys():
                yield from entry_gen(nested[key], _index+(key,), preprocess)

    preproc = kwargs.get('preprocess', None)

    array_labels = kwargs.get('array_labels', {})
    array_params: ArrayParams = {
        "order": [k for k in data_order if k in array_labels],
        "labels": {
            k: parse_ranges(v)
            for k, v in array_labels.items()
        }
    }

    dict_order = [k for k in data_order if k not in array_labels]

    if isinstance(data, np.ndarray):
        return ndarray_to_frame(
            preproc(data) if preproc else data,
            array_params)
    else:
        assert isinstance(data, t.Dict)
        assert len(array_labels) > 0

        indices, data_to_concat = zip(*(
            (index,array)
            for index, array in entry_gen(data, preprocess=preproc)
        ))
        data_to_concat = [ndarray_to_frame(x,array_params) for x in data_to_concat]

        for index, frame in zip(indices,data_to_concat):
            frame[dict_order] = list(index)

        df = pd.concat(data_to_concat)

        df.set_index(dict_order, append=True, inplace=True)

    return df


def load_file(filename: str, **kwargs):

    if filename.endswith(".p"):

        with open(filename, "rb") as f:
            data = pickle.load(f)

        if isinstance(data, pd.DataFrame):
            return data
        elif isinstance(data, t.Dict):
            return dict_to_df(data, **kwargs)
        else:
            raise ValueError(
                (f"Contents of {filename} is of type {type(data)}."
                 "Expecting dictionary or pandas DataFrame."))

    elif filename.endswith(".h5"):
        try:
            return pd.read_hdf(filename)
        except ValueError:
            file = h5py.File(filename)
            return h5_to_df(file, **kwargs)
        raise NotImplementedError("hdf5 not yet implemented")
    else:
        raise ValueError("File must have extension '.p' or '.h5'")


def load_input(filestem: str, replacements: t.Optional[t.Dict] = None,
               regex: t.Optional[t.Dict] = None, **kwargs):
    params: t.Dict = deep_copy_dict(kwargs)

    repl_keys: t.List[str] = formatkeys(filestem)

    str_repl: t.Dict = replacements if replacements else {}
    regex_repl: t.Dict = regex if regex else {}

    assert len(repl_keys) == len(str_repl) + len(regex_repl)
    assert all((
        (k in str_repl or k in regex_repl)
        for k in repl_keys
    ))

    df = []

    actions: t.Dict = params.pop('process', {})
    proc = functools.partial(processor.execute, actions=actions)

    for str_reps, repl_filename in string_replacement_gen(
            filestem, str_repl):
        for reg_reps, regex_filename in file_regex_gen(
                repl_filename, regex_repl):
            logging.debug(f"Loading file: {regex_filename}")
            new_data: pd.DataFrame = proc(load_file(
                filename=regex_filename,
                **params))

            if len(str_reps) != 0:
                new_data[list(str_reps.keys())] = tuple(str_reps.values())
            if len(reg_reps) != 0:
                new_data[list(reg_reps.keys())] = tuple(reg_reps.values())

            df.append(new_data)

    df = pd.concat(df)

    if len(repl_keys) != 0:
        df.set_index(repl_keys, append=True, inplace=True)

    return df


def write_frame(df: pd.DataFrame, filename: str,
                col_to_repl: t.Optional[t.List[str]] = None):
    """Write DataFrame to `filename`. Searches for elements of `col_to_repl`
    in the DataFrame index. Uses first entry of DataFrame for
    keyword replacement in `filename`
    """
    if repl_keys := formatkeys(filename):
        assert len(df) != 0
        assert isinstance(col_to_repl, t.List)
        assert all([k in col_to_repl for k in repl_keys])
        assert all([k in df.index.names for k in repl_keys])

        indices = [
            df.index.names.index(k)
            for k in repl_keys
        ]

        repl = {
            df.index.names[i]: df.iloc[0].name[i]
            for i in indices
        }
        df.reset_index(
            level=indices, drop=True
        ).to_hdf(filename.format(**repl), key='corr')

    else:
        df.to_hdf(filename, key='corr')
