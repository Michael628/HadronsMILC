import os
import numpy as np
import itertools
from string import Formatter
import h5py
import pickle
from dataclasses import dataclass, field
import re
import functools


@dataclass
class FilestemFormatParser:
    """Utility class to parse file strings with multiple keyword replacements

    Attributes
    ----------
    filestem : str
        The file string that will be searched for keyword replacements
    params : dict
        The dictionary from which to find the possible replacements
        for each keyword
    regex : dict, optional
        Provides regex patterns to use as keyword replacements.
        These patterns are used by `traverse_files` to search
        through the file system
    keysort : callable, optional
        A function used to sort the order in which keywords are iterated
        by `traverse_replacements()`
    _current_file : str
        For generator methods `traverse_files`, and `traverse_replacements`
        stores the current file string
    """
    filestem: str
    params: dict = field(default_factory=dict)
    regex: dict = None
    keysort: callable = None

    def __post_init__(self):
        self._current_file = None

        if self.keysort is None:
            self.keysort = self.filestem.index

        self._filekeys = formatkeys(
            self.filestem, key_sort=self.filestem.index)

        self.regex_keyorder = [
            k
            for k in self._filekeys
            if k in self.regex
        ]

        self._filekeys = formatkeys(self.filestem, key_sort=self.keysort)

    @property
    def keys(self):
        return self._filekeys

    def format(self, **kwargs) -> str | functools.partial:
        try:
            return self.filestem.format(**kwargs)
        except KeyError:
            return functools.partial(self.filestem.format, **kwargs)

    def traverse_replacements(self) -> (dict, str):
        """Generator for keyword replacements of `filestem`

        Yields
        ------
        (replacements, filename): (dict, str)
            `replacements` : dict
                The replacement dictionary
                which can be passed to str.format() as kwargs
            `filename` : str
                The `filestem` string formatted by `replacements`
        """

        keys, val_iter = format_iter(
            self.filestem,
            self.params,
            self.keysort
        )

        for vals in val_iter:
            repl = dict(zip(keys, vals))
            self._current_file = self.format(**repl)

            yield repl, self._current_file

        self._current_file = None

    def traverse_regex(self) -> (dict, str):
        """Generator for files found matching regex replacements
        of `_current_file`

        Yields
        ------
        regex_repl : dict
            A dictionary of the matching regex patterns that were matched
        filename : str
            A filename that matches the `_current_file` string using
            the patterns in `regex`, if present
        """

        if self._current_file is None:
            self._current_file = self.filestem

        if len(self.regex) == 0:
            yield {}, self._current_file
        else:
            # Build regex objects to catch each replacement
            regex_repl = {
                k: f"({val})"
                for k, val in self.regex.items()
            }
            try:
                file_pattern = self._current_file.format(**regex_repl)
            except AttributeError:
                file_pattern = self._current_file(**regex_repl)

            # FIX ME: Assumes all regex matches occur in file name,
            # not in the directory path.
            directory, match = os.path.split(file_pattern)

            files: list[str] = os.listdir(directory)

            regex_pattern = re.compile(match)

            for file in files:
                # Skip files that do not match

                regex_matches = regex_pattern.match(file)

                if regex_matches is not None:
                    regex_repl: dict = {
                        key: regex_matches[
                            self.regex_keyorder.index(key)+1
                        ]
                        for key in self.keys
                        if key in self.regex
                    }
                    yield regex_repl, f"{directory}/{file}"


def extractdata(file: str, data_keys: list[str] = None):

    if file[-1] == "p":
        return pickle.load(open(file, "rb"))
    elif file[-2:] == "h5":
        if data_keys is None:
            raise Exception("Name string expected for hdf5 files.")

        with h5py.File(file, "r") as f:
            return np.array(
                [np.array(f[key][()].view(np.complex128), dtype=np.complex128)
                    for key in data_keys]
            )
    else:
        raise Exception("Only pickle and hdf5 files are supported.")


def dictval_iter(input_dict: dict, key_sort=None):
    """Traverses dictionary values for lists of strings and
    builds a list of all combinations of list values

    Parameters
    ----------
    input_dict : dict
        A dictionary who's values can be typed to str or list[str]

    key_sort: Callable, optional
        A function that, if set, is used to sort the key order
    Returns
    -------
    keys, prod: (list[str], iter)
        `keys` is the list of keys of `input_dict`
        `prod` is a list of all combinations of the values of `input_dict`
    """

    keys = list(input_dict.keys())
    if key_sort:
        keys.sort(key=key_sort)

    prod_list = [
                    [str(item) for item in input_dict[k]]
                    if type(input_dict[k]) is list
                    else [str(input_dict[k])]
                    for k in keys
                ]

    return keys, itertools.product(*prod_list)


def getvalue(dict_in: dict, path: list[str] | str) -> dict:

    _path = [path] if isinstance(path, str) else path

    val = dict_in
    try:
        for key in _path:
            val = val[key]
    except KeyError:
        return None
    return val


def findvalues(dict_in: dict, paths: list[str], keys: list[str]):
    """Searches the `paths` of `dict_in` for `keys and returns
    a dictionary of the values found for each key. Returns the
    first value found when traversing `paths`
    """

    dict_out = {}
    for path in paths:
        for key, val in zip(keys, map(
                    lambda x: getvalue(dict_in, x),
                    ([path]+[key] for key in keys))):

            if val is not None:
                dict_out[key] = val

    if any(k not in dict_out for k in keys):
        raise KeyError(
            " ".join((
                "Could not find values for key(s):",
                ", ".join((
                    k for k in keys
                    if k not in dict_out
                )))
            )
        )

    return dict_out


def formatkeys(format_string: str, key_sort: callable = None) -> list[str]:
    key_list = list({
        k[1]
        for k in Formatter().parse(format_string)
        if k[1] is not None
    })

    if key_sort is not None:
        key_list.sort(key=key_sort)

    return key_list


def formatdict(format_string: str, source_dict: dict) -> dict:
    """Builds a dictionary from elements of `source_dict` of required variables
    in `format_string`.

    Parameters
    ----------
    format_string : str
        A string with named format variables that need replacement

    source_dict : dict
        The dictionary from which the returned dict is built.

    Returns
    -------
    sub_dict : dict
        Dictionary entries of `source_dict` matching all keyword elements found
        in `format_string`.
    """

    key_list = formatkeys(format_string)

    sub_dict = {}
    for key in key_list:
        if key in source_dict:
            sub_dict[key] = source_dict[key]

    return sub_dict


def format_iter(format_string: str, source_dict: dict,
                key_sort=None) -> list[tuple]:
    """A list of all combinations of format parameters designated in
    `format_string` and found in `source_dict`
    """

    return dictval_iter(formatdict(format_string, source_dict), key_sort)


def setdictval(dict_out: dict, keys: list[str],
               value=field(default_factory=dict), overwrite=False):
    """Recursively builds nested dictionary. Nesting is in order of `keys`
    """

    if len(keys) == 1:
        if overwrite or keys[0] not in dict_out:
            dict_out[keys[0]] = value
    else:
        inner_dict = dict_out.get(keys[0], {})
        setdictval(inner_dict, keys[1:], value)
        dict_out.update({keys[0]: inner_dict})


def dict_to_corr(corr_dict):
    """Traverses a nested dictionary and returns a multidimensional numpy array
    """

    def dict_to_corr_helper(index, corr_out, corr_dict):
        # If we've reached the lowest level of the dictionary,
        # copy contents and increment counter
        if type(corr_dict) is not dict:
            corr_out[index, :] = corr_dict
            index += 1
            return index

        # Otherwise recurse over next lower level in dictionary
        for k, v in sorted(corr_dict.items()):
            print(f"time: {k}")
            index = dict_to_corr_helper(index, corr_out, corr_dict[k])
        return index

    dict_shape = tuple()

    # Determine shape of correlator based on number of keys in
    # each nested dictionary
    d = corr_dict
    while type(d) is dict:
        dict_shape = dict_shape + (len(d.keys()),)
        d = next(iter(d.values()))
    inner_shape = d.shape

    dict_size = np.prod(dict_shape)
    corr_out = np.empty((dict_size,)+inner_shape, dtype=np.complex128)

    # Recursively fill correlator with nested dictionary entries
    index = 0
    for v in corr_dict.values():
        index = dict_to_corr_helper(index, corr_out, v)

    # Reshape to multidimensional structure
    corr_out = np.reshape(corr_out, dict_shape+inner_shape)
    return corr_out
