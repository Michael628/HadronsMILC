import yaml
import copy
from collections.abc import Mapping


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


def load_params_join(YAMLEns, YAMLAll):
    """Concatenate two YAML parameter files and load
    We need this because YAMLEns defines a reference needed
    by YAMLAll"""

    # Initial parameter file
    # try:
    ens = open(YAMLEns, 'r').readlines()
    all = open(YAMLAll, 'r').readlines()
    param = yaml.safe_load("".join(ens + all))
    # except:
    # print("ERROR: Error loading the parameter files", YAMLEns, YAMLAll)
    # sys.exit(1)

    return param


############################################################
def update_param(param, param_update):
    """Update the param dictionary according to terms in param_update"""

    # Updating is recursive in the tree so we can update selected branches
    # leaving the remainder untouched
    for b in param_update.keys():
        try:
            k = param_update[b].keys()
            n = len(k)
        except AttributeError:
            n = 0

        if b in param.keys() and n > 0:
            # Keep descending until we run out of branches
            update_param(param[b], param_update[b])
        else:
            # Then stop, replacing just the last branch or creating a new one
            param[b] = param_update[b]

    return param
