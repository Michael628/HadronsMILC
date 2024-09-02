import numpy as np
import itertools
from string import Formatter
import h5py
import pickle
from dataclasses import field


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
    keys, prod: (list[str], list[tuple[str]])
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
        inner_dict = {} if keys[0] not in dict_out else dict_out[keys[0]]
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
            index = dict_to_corr_helper(index, corr_out, corr_dict[k])
        return index

    dict_shape = tuple()

    # Determine shape of correlator based on number of keys in
    # each nested dictionary
    d = corr_dict
    while type(d) is dict:
        dict_shape = dict_shape + (len(d.keys()),)
        d = list(d.values())[0]
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
