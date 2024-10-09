import numpy as np
import h5py
import pickle
from dataclasses import dataclass, field
import logging
import copy
import itertools
import pandas as pd
from time import perf_counter
from python_scripts.processing.format \
    import FilestemFormatParser as FormatParser


DF_INTEGER_KEYS = ['cfg', 'time']


@dataclass
class DataLoad:
    """Loads data from hdf5 or pickled dictionaries of ndarrays
    and collects them into a pandas data frame

    Attributes
    ----------
    input_params : dict
        kwargs to forward to Fo
    data_keys : list[str]
        The keys
    data_labels : dict, optional
        Required for pickled dictionary files. See `extractdata` docstring.
    datapaths : dict, optional
        Required for hdf5 files. See `extractdata` docstring.
    *overwrite : bool = False
        When false and an existing output file is found, skips input files
        that correspond to data already in existing output file
    """
    input_params: dict
    data_keys: list[str]
    data_labels: dict = field(default_factory=dict)
    datapaths: dict = field(default_factory=dict)
    preprocess: callable = None

    @classmethod
    def traverse_dict(cls, result, source: dict | np.ndarray, key_labels: list,
                      depth: int = 0):
        """Generator, which traverses `source`, a nested dictionary, and
        yields a dictionary keyed by `key_labels` with values corresponding
        to the key found at each depth level in `source` (in the order
        they appear in `key_labels`). The 'data' key in the returned dictionary
        holds the ndarray found at max depth

        Parameters
        ----------
        result : dict
            dictionary containing `key_labels` + ['data']. The value
            corresponding to 'data' is the ndarray object found at max depth
            of the nested dictionary
        """
        label = key_labels[depth]
        if label not in result:
            result[label] = []
        if 'data' not in result:
            result['data'] = []

        if not isinstance(next(iter(source.values())), dict):
            result[label] += list(source.keys())
            result['data'] += list(source.values())
        else:
            for key, val in source.items():
                cls.traverse_dict(result, val, key_labels, depth+1)
                newdata_len = len(result['data']) - len(result[label])
                result[label] += [key]*newdata_len

    @classmethod
    def extractdata(cls, filename: str,  key_labels: list[str],
                    datapaths: dict = None, data_labels: dict = None,
                    preprocess=None):
        """Pull data from file (pickle or hdf5)

        Parameters
        ----------
        filename : str
            The file that will be opened.
        key_labels: list[str]
            Labels for data elements found in `file`
        data_labels: dict, optional
            Required for pickled dictionaries. Used to flatten multidimensional
            ndarrays
        datapaths: dict, optional
            Required for hdf5 files.

        Yields
        -------
        total : dict
            A dictionary that conforms
            to the information pulled from the file.

        Notes
        -----
        *pickle*
        If `file` is a pickled object file, then `data_keys` is a list[str]
        where the entries of `data_keys` are labels for a path of keys
        through nested dictionaries

        *hdf5*
        If `file` is of type hdf5, then `data_keys` is a dict
        where entries in `data_keys` have the structure:
        'object_path':{
            'key': 'h5_label'
        }
        Thus, `extractdata` returns a list of dicts with items:
            - (key, f['object_path']['h5_label']) or
            - (key, f['object_path'].attrs['h5_label'])
        """
        if filename.endswith(".p"):
            file = pickle.load(open(filename, "rb"))

            result = {}
            cls.traverse_dict(result, file, key_labels)

            inner_keys = [k for k in key_labels if k not in result]
            inner_vals = [data_labels[k]
                          for k in key_labels if (k in data_labels)]

            # Convert 'num1..num2' to integer range
            inner_vals = [
                range(*map(int, v.split("..")))
                if isinstance(v, str)
                else v
                for v in inner_vals
            ]

            file_data = result.pop('data')

            if not all([len(inner_vals[i]) == file_data[0].shape[i]
                        for i in range(len(inner_vals))]):
                raise ValueError(
                    "data_labels dimensions do not match file data.")

            if preprocess is not None:
                file_data = map(preprocess, file_data)

            file_data = np.concatenate([x[np.newaxis] for x in file_data])
            data_shape = file_data.shape[len(inner_vals)+1:]
            file_data = file_data.reshape((-1,)+data_shape)

            # Postrocessing keys
            if 'series.cfg' in result:
                series, cfg = zip(*[
                    x.split('.')
                    for x in result.pop('series.cfg')
                ])
                result['series'] = series
                result['cfg'] = cfg
            if 'time' in result and result['time'][0].startswith('t'):
                result['time'] = [x[1:] for x in result['time']]

            indices = itertools.product(
                list(zip(*result.values())),
                itertools.product(*inner_vals)
            )
            indices = [sum(x, ()) for x in indices]
            columns = list(result.keys())+inner_keys+['data']

            items = (
                index + (file_data[i],)
                for i, index in enumerate(indices)
            )
            total = pd.DataFrame(items, columns=columns)
            yield total

        elif filename.endswith(".h5"):
            assert datapaths is not None
            result = []
            with h5py.File(filename, "r") as file_data:
                for i, (dpath, ditem) in enumerate(datapaths.items()):
                    result.append({})
                    for key, val in ditem.items():
                        try:  # Try to interpret as attribute
                            result[i][key] = \
                                file_data[dpath].attrs[val][0].decode('utf-8')
                        except KeyError:  # Try to interpret as ndarray
                            result[i][key] = np.array(
                                file_data[dpath][val][()].view(np.complex128),
                                dtype=np.complex128)
            for item in result:
                yield item
        else:
            raise ValueError("File must have extension '.p' or '.h5'")

    def __post_init__(self):

        self._df = None

        self._input_parser = FormatParser(**self.input_params)

        full_keyset: list = self._input_parser.keys + self.data_keys

    def execute(self) -> dict:

        res = {}

        for i, (file_repl, infile) in \
                enumerate(self._input_parser.traverse_files()):

            logging.debug(f"Processing file: {infile}")

            for new_data in self.extractdata(infile,
                                             datapaths=self.datapaths,
                                             key_labels=self.data_keys,
                                             data_labels=self.data_labels,
                                             preprocess=self.preprocess):

                if self._df is None:
                    self._df = pd.DataFrame(columns=list(
                        file_repl.keys())+list(new_data.keys()))

                self._df = pd.concat([self._df, new_data], ignore_index=True)
                self._df.fillna(file_repl, inplace=True)

        intkeys = [k for k in self._df.columns if k in DF_INTEGER_KEYS]
        self._df[intkeys] = self._df[intkeys].apply(pd.to_numeric)

        # if len(data_index) != 0:
        # self._df.set_index(data_index, inplace=True)

        return self._df
