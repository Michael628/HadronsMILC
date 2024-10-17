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
