import concurrent.futures
import copy
import itertools
import logging
import os
import re
import typing as t
from collections.abc import Mapping
from functools import partial
from string import Formatter

import yaml

import python_scripts as ps

procFn = t.Callable[[str, t.Any], t.Any]


class ReadOnlyDict(Mapping):
    def __init__(self, initial_data):
        self._data = dict(initial_data)

    def __getitem__(self, key):
        return self._data[key]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def __repr__(self):
        return f"{self.__class__.__name__}({self._data})"

    def __setitem__(self, key, value):
        raise TypeError("This dictionary is read-only")

    def __delitem__(self, key):
        raise TypeError("This dictionary is read-only")


def deep_copy_dict(read_dict):
    def is_dict(d):
        return isinstance(d, ReadOnlyDict) or \
            isinstance(d, dict)

    """Recursively create a writable deep copy of dict."""
    if not is_dict(read_dict):
        raise ValueError("Input must be a ReadOnlyDict or dict instance")

    writable_copy = {}
    for key, value in read_dict.items():
        if is_dict(value):
            # Recursively copy nested ReadOnlyDict
            writable_copy[key] = deep_copy_dict(value)
        else:
            # For other types, just copy the value
            writable_copy[key] = copy.deepcopy(value)

    return writable_copy


def load_param(file):
    """Read the YAML parameter file"""

#    try:
    param = yaml.safe_load(open(file, 'r'))
#    except subprocess.CalledProcessError as e:
#        print("WARNING: load_param failed for", e.cmd)
#        print("return code", e.returncode)
#        sys.exit(1)

    return param


# ------ Format/File Search Functions ------ #
def process_params(**params) -> t.Dict:
    """Does nothing to `val` parameters of type list, but converts
    strings to number ranges.
    """

    param_out = deep_copy_dict(params)

    for key, val in params.items():

        if isinstance(val, str):
            if '..' in val:
                range_input = list(map(int, val.split("..")))
                range_input[1] += 1
                param_out[key] = list(range(*range_input))
            else:
                param_out[key] = [val]
        elif not isinstance(val, t.List):
            logging.debug(f'Removing key: {key} from parameters.')
            param_out.pop(key)

    return param_out


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


def process_files(filestem: str, processor: procFn,
                  replacements: t.Optional[t.Dict] = None,
                  regex: t.Optional[t.Dict] = None) -> t.List:

    repl_keys: t.List[str] = formatkeys(filestem)

    str_repl: t.Dict = replacements if replacements else {}
    regex_repl: t.Dict = regex if regex else {}

    logging.debug(f'repl_keys: {sorted(repl_keys)}')
    logging.debug(f'str_repl keys: {sorted(str_repl.keys())}')
    logging.debug(f'regex_repl keys: {sorted(regex_repl.keys())}')
    assert len(repl_keys) == len(str_repl) + len(regex_repl)
    assert all(((k in str_repl or k in regex_repl) for k in repl_keys))

    collection: t.List = []

    def file_gen():
        for str_reps, repl_filename in string_replacement_gen(
                filestem, str_repl):
            for reg_reps, regex_filename in file_regex_gen(
                    repl_filename, regex_repl):

                str_reps.update(reg_reps)
                yield regex_filename, deep_copy_dict(str_reps)

    for filename, reps in file_gen():
        try:
            new_result = processor(filename, reps)
        except StopIteration as e:
            if e.args:
                assert len(e.args) == 1
                collection.append(e.args[0])
            break
        collection.append(new_result)

    return collection
# ------ End Format/File Search Functions ------ #


