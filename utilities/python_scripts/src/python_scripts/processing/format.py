import os
import itertools
from string import Formatter
from dataclasses import dataclass, field
import re
import functools
import typing as t


@dataclass
class FilestemFormatBase:
    @staticmethod
    def getvalue(haystack: t.Dict, needle: t.List[str]) -> t.Dict:
        """Traverse `needle` in nested dictionary `haystack`

        Returns
        -------
        res : Any
            The value that is found or None if the path contains
            an invalid key. Returns `haystack` if `needle` is empty.
        """
        path = [needle] if isinstance(needle, str) else needle

        res = haystack
        try:
            for key in path:
                res = res[key]
        except KeyError:
            res = None

        return res

    @staticmethod
    def _fix_input(s):
        """Turn input into type list[list[str]]"""

        res = s
        if isinstance(s, str):
            res = [[s]]
        elif isinstance(s, t.List):
            if len(s) == 0 or isinstance(s[0], str):
                res = [s]

        if len(res[0]) != 0 and not isinstance(res[0][0], str):
            raise ValueError("Expecting str, list[str], or list[list[str]].")

        return res

    @classmethod
    def getvalues(cls, haystack: t.Dict,
                  paths: t.List[t.List[str]], needles: t.List[t.List[str]]):
        """Searches the `paths` of `haystack` for each key in `needles`.
        Takes the first value found when traversing `paths`

        Returns
        -------
        res : dict
            A dictionary of the values found for each key in `needle`
        """

        ps = cls._fix_input(paths)
        ns = cls._fix_input(needles)

        res = {}
        for needle in ns:
            for p, substack in zip(
                ps, map(lambda x: cls.getvalue(haystack, x), ps)
            ):
                try:
                    if isinstance(substack, dict):
                        res[needle] = substack[needle]
                        break
                except KeyError:
                    pass

            if needle not in res:
                raise KeyError(f"Could not find value for key '{needle}'")

        return res

    @staticmethod
    def formatkeys(format_string: str, keysort: t.Callable = None) -> t.List[str]:
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

    @classmethod
    def formatdict(cls, format_string: str, **kwargs) -> t.Dict:
        """Builds a dictionary from elements of `source_dict` of required
        variables in `format_string`.

        Parameters
        ----------
        format_string : str
            A string with named format variables that need replacement

        source_dict : dict
            The dictionary from which the returned dict is built.

        Returns
        -------
        sub_dict : dict
            Dictionary entries of `source_dict` matching all keyword elements
            found in `format_string`.
        """

        key_list = cls.formatkeys(format_string)

        sub_dict = {}
        for key in key_list:
            if key in kwargs:
                sub_dict[key] = kwargs[key]

        return sub_dict

    @staticmethod
    def dictval_iter(input_dict: t.Dict, keysort=None) -> (t.List[str], iter):
        """Traverses dictionary values for lists of strings and
        builds a list of all combinations of list values

        Parameters
        ----------
        input_dict : dict
            A dictionary who's values can be typed to str or list[str]

        keysort: Callable, optional
            A function that, if set, is used to sort the key order
        Returns
        -------
        keys, prod: (list[str], iter)
            `keys` is the list of keys of `input_dict`
            `prod` is a list of all combinations of the values of `input_dict`
        """

        keys = list(input_dict.keys())
        if keysort:
            keys.sort(key=keysort)

        prod_list = [
            [str(item) for item in input_dict[k]]
            if isinstance(input_dict[k], t.List)
            else [str(input_dict[k])]
            for k in keys
        ]

        return keys, itertools.product(*prod_list)

    @classmethod
    def format_iter(
        cls, format_string: str, source_dict: t.Dict, keysort=None
    ) -> t.List[t.Tuple]:
        """A list of all combinations of format parameters designated in
        `format_string` and found in `source_dict`
        """

        return cls.dictval_iter(
            cls.formatdict(format_string, **source_dict), keysort
        )

    @staticmethod
    def setdictval(
        dict_out: t.Dict,
        keys: t.List[str],
        value=None,
        overwrite=False,
    ):
        """Recursively builds nested dictionary.
        Nesting is in order of `keys`
        """

        if value is None:
            value = {}

        if len(keys) == 1:
            if overwrite or keys[0] not in dict_out:
                dict_out[keys[0]] = value
        else:
            inner_dict = dict_out.get(keys[0], {})
            setdictval(inner_dict, keys[1:], value)
            dict_out.update({keys[0]: inner_dict})


@dataclass
class FilestemFormatParser(FilestemFormatBase):
    """Utility class to parse file strings with multiple keyword replacements

    Attributes
    ----------
    filestem : str
        The file string that will be searched for keyword replacements
    replacements : dict
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
    replacements: t.Dict = field(default_factory=dict)
    regex: t.Dict = field(default_factory=dict)
    keysort: t.Callable = None

    def __post_init__(self):
        self._current_file = None

        if self.keysort is None:
            self.keysort = self.filestem.index
        self._filekeys = self.formatkeys(self.filestem, keysort=self.keysort)

        filestem_keyorder = self.formatkeys(
            self.filestem, keysort=self.filestem.index
        )
        self.regex_keyorder = [k for k in filestem_keyorder if k in self.regex]

    @property
    def keys(self):
        return self._filekeys

    def format(self, **kwargs) -> t.Union[str, functools.partial]:
        try:
            return self.filestem.format(**kwargs)
        except KeyError:
            return functools.partial(self.filestem.format, **kwargs)

    def traverse_replacements(self) -> (t.Dict, str):
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

        keys, val_iter = self.format_iter(
            self.filestem, self.replacements, self.keysort
        )

        for vals in val_iter:
            repl = dict(zip(keys, vals))
            self._current_file = self.format(**repl)

            yield repl, self._current_file

        self._current_file = None

    def traverse_regex(self) -> (t.Dict, str):
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
            regex_repl = {k: f"({val})" for k, val in self.regex.items()}
            try:
                file_pattern = self._current_file.format(**regex_repl)
            except AttributeError:
                file_pattern = self._current_file(**regex_repl)

            # FIX ME: Assumes all regex matches occur in file name,
            # not in the directory path.
            directory, match = os.path.split(file_pattern)

            files: t.List[str] = os.listdir(directory)

            regex_pattern = re.compile(match)

            for file in files:
                # Skip files that do not match

                regex_matches = regex_pattern.match(file)

                if regex_matches is not None:
                    regex_repl: dict = {
                        key: regex_matches[self.regex_keyorder.index(key) + 1]
                        for key in self.keys
                        if key in self.regex
                    }
                    yield regex_repl, f"{directory}/{file}"

    def traverse_files(self) -> (t.Dict, str):
        full_repl = {}
        for replacement, _ in self.traverse_replacements():
            full_repl.update(replacement)
            for regex_repl, filename in self.traverse_regex():
                full_repl.update(regex_repl)
                yield full_repl, filename
