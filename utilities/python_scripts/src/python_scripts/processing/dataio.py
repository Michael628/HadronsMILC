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
from python_scripts.nanny.todo_utils import load_param
from python_scripts.utils import deep_copy_dict
DF_INTEGER_KEYS = ['cfg', 'time', 'dt']


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
        yield {}, functools.partial(lambda: fstring)
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
        array_order: t.List[str] = ['idx'],
        array_labels: t.Dict = {}) -> pd.DataFrame:
    """Converts ndarray into a pandas.Series object indexed by the values of
    `array_labels`, with indices named by the keys of `array_labels.
    `label_order` provides the ordering of the labels in the multidimensional
    ndarray.
    """

    if len(array_order) == 1 and array_order[0] == 'idx':
        array_labels['idx'] = list(range(np.prod(array.shape)))

    assert len(array_labels) == len(array_order)

    indices = [array_labels[k] for k in array_order]

    index: pd.MultiIndex = pd.MultiIndex.from_tuples(
        itertools.product(*indices),
        names=array_order
    )

    return pd.Series(
        array.reshape((-1,)),
        index=index,
        name='corr'
    ).to_frame()


def load_pickle(
        filename: str,
        data_order: t.List[str] = [],
        **kwargs) -> t.Any:

    def do_nothing(x):
        return x

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

    def entry_gen(
        nested: t.Dict,
        _index: t.List = [],
        **kwargs) \
            -> t.Generator[t.Tuple[t.Tuple, pd.DataFrame], None, None]:
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
                    tuple(_index+[key]),
                    ndarray_to_frame(preproc(val), **kwargs)
                )
        else:
            for key in nested.keys():

                _index.append(key)
                yield from entry_gen(nested[key], _index, **kwargs)
                _index.pop()

    with open(filename, "rb") as f:
        data = pickle.load(f)

    preproc = kwargs.get('preprocess', do_nothing)

    array_labels = kwargs.get('array_labels', {})
    array_labels = {
        k: parse_ranges(v)
        for k, v in array_labels.items()
    }

    array_order = [k for k in data_order if k in array_labels]
    dict_order = [k for k in data_order if k not in array_labels]

    if isinstance(data, np.ndarray):
        return ndarray_to_frame(preproc(data), array_order, array_labels)
    else:
        assert isinstance(data, t.Dict)
        assert len(array_labels) > 0

        df = None
        for outer, frame in entry_gen(
                data,
                array_order=array_order,
                array_labels=array_labels):
            for label, key in zip(dict_order, outer):
                frame[label] = key
            if df is None:
                df = frame
            else:
                df = pd.concat([df, frame])

        assert isinstance(df, pd.DataFrame)

        df.set_index(dict_order, append=True, inplace=True)

    return df


def load_file(filename: str, **kwargs):

    if filename.endswith(".p"):
        return load_pickle(filename, **kwargs)
    elif filename.endswith(".h5"):
        raise NotImplementedError("hdf5 not yet implemented")
    else:
        raise ValueError("File must have extension '.p' or '.h5'")


def load_input(config: t.Dict, **kwargs):
    params: t.Dict = deep_copy_dict(kwargs)

    repl_keys: t.List[str] = formatkeys(config['filestem'])
    str_repl: t.Dict = config.get("replacements", {})

    regex_repl: t.Dict = config.get('regex', {})

    assert len(repl_keys) == len(str_repl) + len(regex_repl)
    assert all((
        (k in str_repl or k in regex_repl)
        for k in repl_keys
    ))

    df = None

    actions: t.Dict = params.pop('process', {})
    proc = functools.partial(processor.execute, actions=actions)

    for str_reps, repl_filename in string_replacement_gen(
            config['filestem'], str_repl):
        for reg_reps, regex_filename in file_regex_gen(
                repl_filename, regex_repl):

            new_data: pd.DataFrame = proc(load_pickle(
                filename=regex_filename,
                **params))

            if len(str_reps) != 0:
                new_data[list(str_reps.keys())] = tuple(str_reps.values())
            if len(reg_reps) != 0:
                new_data[list(reg_reps.keys())] = tuple(reg_reps.values())

            if df is None:
                df = new_data
            else:
                df = pd.concat([df, new_data])

    assert isinstance(df, pd.DataFrame)

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
