import asyncio
from concurrent.futures import ThreadPoolExecutor

import python_scripts as ps

import os
import logging
import itertools
from functools import partial
import numpy as np
import pandas as pd
import h5py
import typing as t

from python_scripts import utils
from python_scripts.processing import (
    processor,
    config
)


dataFrameFn = t.Callable[[np.ndarray], pd.DataFrame]
loadFn = t.Callable[[str], pd.DataFrame]


# ------ Data structure Functions ------ #
def ndarray_to_frame(
        array: np.ndarray,
        array_params: config.LoadArrayConfig) -> pd.DataFrame:
    """
    Converts a multidimensional numpy array into a pandas DataFrame with a MultiIndex.

    Parameters:
    array : np.ndarray
        The input numpy array to be transformed into a DataFrame.
    array_params : config.LoadArrayConfig
        The configuration object that provides the order of dimensions (order),
        and labels for indexing each dimension.

    Returns:
    pd.DataFrame
        A DataFrame where the rows are indexed by a MultiIndex created based on
        the combinations of labels defined in array_params, and the single column
        'corr' holds the values from the input array.

    Raises:
    AssertionError
        If the number of label sets in array_params does not match the number
        of dimensions defined in array_params.order.

    Behavior:
    - If the array_params.order contains only 'dt', the function generates
      sequential labels for the 'dt' dimension.
    - Constructs a MultiIndex using the Cartesian product of the label sets
      specified in array_params.labels, adhering to the order in array_params.order.
    - Flattens the input array and pairs its values to the MultiIndex, creating
      a single-column DataFrame.
    """

    if len(array_params.order) == 1 and array_params.order[0] == 'dt':
        array_params.labels['dt'] = list(range(np.prod(array.shape)))

    assert len(array_params.labels) == len(array_params.order)

    indices = [array_params.labels[k] for k in array_params.order]

    index: pd.MultiIndex = pd.MultiIndex.from_tuples(
        itertools.product(*indices),
        names=array_params.order
    )

    return pd.Series(
        array.reshape((-1,)),
        index=index,
        name='corr'
    ).to_frame()


def h5_to_frame(file: h5py.File,
                data_to_frame: t.Dict[str, dataFrameFn],
                h5_params: config.LoadH5Config) -> pd.DataFrame:
    """
    Converts datasets from an HDF5 file to a Pandas DataFrame.

    Parameters:
    file : h5py.File
        The HDF5 file object containing the datasets to be converted.
    data_to_frame : Dict[str, Callable]
        A dictionary mapping dataset labels to functions that convert their data into DataFrames.
    h5_params : config.LoadH5Config
        Configuration object specifying parameters for loading datasets from the HDF5 file,
        including dataset names and additional metadata.

    Returns:
    pd.DataFrame
        A Pandas DataFrame constructed from the datasets in the HDF5 file. If more than one dataset
        is included, an additional index column is added to distinguish between them.

    Raises:
    ValueError
        If a specified dataset is not found in the HDF5 file.
    """
    assert all(k in data_to_frame.keys() for k in h5_params.datasets.keys())

    df = []
    for k, v in h5_params.datasets.items():
        
        if isinstance(v,str):
            dataset_label = iter((v,) if v in file else ())
        else:
            dataset_label = (x for x in v if x in file)
        
        try:
            data = file[next(dataset_label)][:].view(np.complex128)
        except StopIteration:
            raise ValueError(f"dataset {k} not found in file.")

        frame = data_to_frame[k](data)
        if len(h5_params.datasets) > 1:
            frame[h5_params.name] = k
        df.append(frame)

    df = pd.concat(df)

    if len(h5_params.datasets) > 1:
        df.set_index(h5_params.name, append=True, inplace=True)

    return df


def frame_to_dict(df: pd.DataFrame, dict_depth: int) \
        -> t.Union[t.Dict, np.ndarray]:
    """
    Converts a pandas DataFrame into a dictionary or a numpy array depending on the specified depth.

    Parameters:
    df : pandas.DataFrame
        Input DataFrame to be converted.
    dict_depth : int
        Depth of the dictionary to create. If 0, a numpy array is returned instead of a dictionary.

    Returns:
    Union[Dict, numpy.ndarray]
        If `dict_depth` is 0, returns a multidimensional numpy array reshaped based on the levels of the index.
        Otherwise, returns a dictionary keyed by the concatenated indices up to the `dict_depth` level, and values
        are reshaped numpy arrays based on the remaining index levels.

    Behavior:
        Ensures that the `dict_depth` is within the permissible range (0 to the number of index levels).
        Reshapes the values of the DataFrame (assumed column name 'corr') into a multidimensional numpy array
        based on the levels of the index exceeding the specified `dict_depth`.
        Concatenates multi-level index keys up to the specified depth into a string format using '.' as a separator
        when returning a dictionary.

    Important Notes:
        - The DataFrame is expected to have a multi-index and a column named 'corr'.
        - If dict_depth exceeds the number of index levels or is negative, the function raises an assertion error.
        - The DataFrame is sorted based on the index before processing to ensure consistency in the output.
    """
    num_indices = len(df.index.names)
    assert dict_depth >= 0
    assert dict_depth <= num_indices

    shape = [
        len(df.index.get_level_values(i).drop_duplicates())
        for i in range(dict_depth, num_indices)
    ]
    shape = tuple([-1] + shape) if dict_depth != 0 else tuple(shape)

    keys = [
        df.sort_index().index.get_level_values(i).drop_duplicates().to_list()
        for i in range(dict_depth)
    ]

    def join_str_fn(x):
        return ".".join(map(str, x))

    keys = list(map(join_str_fn, list(itertools.product(*keys))))

    array = df.sort_index()['corr'].to_numpy().reshape(shape)

    if dict_depth == 0:
        return array
    else:
        return {k: array[i] for k, i in zip(keys, range(len(array)))}


def dict_to_frame(
        data: t.Union[t.Dict, np.ndarray],
        data_to_frame: dataFrameFn,
        dict_labels: t.Tuple[str] = ()) -> pd.DataFrame:
    """
    Processes nested dictionaries or NumPy arrays and converts them into a pandas DataFrame.

    Parameters:
    data: Dict or ndarray
        The main input which could either be a nested dictionary structure containing NumPy arrays
        as values or a NumPy array itself.
    data_to_frame: callable
        A function used to convert an input ndarray into a pandas DataFrame.
    dict_labels: tuple, optional
        A tuple containing strings which represent labels for additional columns created based on
        dictionary keys in the nested structure.

    Returns:
    DataFrame
        A pandas DataFrame created from the input data. If `data` is a dictionary, its indices
        (keys in hierarchical order) are appended as additional columns using the strings from
        `dict_labels`. If `data` is an ndarray, the function processes it directly via `data_to_frame`.

    Raises:
    AssertionError
        If `data` is of an unsupported type, i.e., not a dictionary or NumPy array, or if values in
        a nested dictionary are inconsistent types.
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
# ------ End data structure functions ------ #


# ------ Input functions ------ #
async def load_files(filestem: str, file_loader: loadFn,
               replacements: t.Optional[t.Dict] = None,
               regex: t.Optional[t.Dict] = None):
    """
    Loads files and processes them using a provided processing function and optional replacements or regex operations.

    Parameters:
    filestem (str): The base name or stem of the files to be loaded and processed.
    file_loader (loadFn): A callable function responsible for file loading.
    replacements (t.Optional[t.Dict]): Optional dictionary containing key replacements to apply to `filestem`.
    regex (t.Optional[t.Dict]): Optional dictionary containing regex patterns as key replacements to `filestem`, matching
     all files found on disk according to the resulting pattern.

    Returns:
    pd.DataFrame: A pandas DataFrame result after processing the loaded files.
    """
    async def proc(filename: str, repl: t.Dict) -> pd.DataFrame:
        logging.debug(f"Loading file: {filename}")
        new_data: pd.DataFrame = file_loader(filename)

        if len(repl) != 0:
            new_data[list(repl.keys())] = tuple(repl.values())

        return new_data

    files = utils.process_files(filestem, lambda *x: x, replacements, regex)

    async_tasks: t.List[asyncio.Task]
    async with asyncio.TaskGroup() as tg:
        async_tasks = [tg.create_task(proc(*file)) for file in files]

    return [t.result() for t in async_tasks]


def load(io_config: config.DataioConfig) -> pd.DataFrame:
    """
    Load data from a variety of file formats and process it into a pandas DataFrame.

    Parameters:
        io_config (config.DataioConfig): Configuration object containing file and processing parameters.

    Returns:
        pd.DataFrame: A DataFrame containing the processed data.

    Raises:
        ValueError: If the file's contents are of an unsupported type or the file extension is unrecognized.

    The function supports loading data from pickle (.p or .npy) and HDF5 (.h5) formatted files.
    For pickle files, it uses a custom loader that converts ndarrays or dictionaries into DataFrames.
    For HDF5 files, it reads the data using pandas' HDF support or processes it manually with custom configurations.
    Actions specified in the io_config object are then applied to the processed data.
    """
    def pickle_loader(filename: str):
        """
        Load data from a given file into a pandas DataFrame.

        Parameters:
            filename (str): File path to the pickle file.

        Returns:
            pd.DataFrame: Data loaded from the file and transformed into a
                          pandas DataFrame.

        Raises:
            ValueError: If the file contents are not a dictionary or pandas DataFrame.

        Notes:
            - The file specified by 'filename' is expected to be a numpy file (.npy)
              with either a dictionary or pandas DataFrame.
            - Data transformation uses the settings provided in 'io_config' to
              ensure compatibility with the DataFrame structure.
        """
        dict_labels: t.Tuple = tuple(io_config.dict_labels)

        array_params: config.LoadArrayConfig = io_config.array_params

        data_to_frame = partial(ndarray_to_frame, array_params=array_params)

        data = np.load(filename, allow_pickle=True)
        if isinstance(data, np.ndarray) and len(data.shape) == 0:
            data = data.item()

        if isinstance(data, pd.DataFrame):
            return data
        elif isinstance(data, np.ndarray):
            return dict_to_frame(data, data_to_frame=data_to_frame)
        elif isinstance(data, t.Dict):
            return dict_to_frame(data,
                                 data_to_frame=data_to_frame,
                                 dict_labels=dict_labels)
        else:
            raise ValueError(
                (f"Contents of {filename} is of type {type(data)}."
                 "Expecting dictionary or pandas DataFrame."))

    def h5_loader(filename: str):
        """
        Load data based on the given I/O configuration.

        This function attempts to load data from an HDF5 file.
        If the file cannot be loaded directly using pandas due to certain errors (e.g., ValueError, NotImplementedError),
        it will fallback to using the provided H5 loading configurations from the `io_config`.
        These configurations specify parameters for handling data arrays and conversion methods for HDF5 structures.

        Parameters:
            filename (str): File path to the pickle file.

        Returns:
            pd.DataFrame
                A DataFrame containing formatted data from the specified HDF5 file.
        """
        try:
            return pd.read_hdf(filename)
        except (ValueError, NotImplementedError):
            assert io_config.h5_params is not None

            h5_params: config.LoadH5Config = io_config.h5_params

            array_params: t.Dict[str, config.LoadArrayConfig]
            array_params = io_config.array_params

            data_to_frame = {
                k: partial(ndarray_to_frame, array_params=array_params[k])
                for k in array_params.keys()
            }
            file = h5py.File(filename)

            return h5_to_frame(file, data_to_frame, h5_params)

    replacements: t.Dict = io_config.replacements
    regex: t.Dict = io_config.regex
    filestem: str = io_config.filestem

    if filestem.endswith(".p") or filestem.endswith(".npy"):
        file_loader = pickle_loader
    elif filestem.endswith(".h5"):
        file_loader = h5_loader
    else:
        raise ValueError("File must have extension '.p' or '.h5'")

    return load_files(
        filestem, file_loader, replacements, regex
    )

# ------ End Input functions ------ #


# ------ Input functions ------ #
def write_data(df: pd.DataFrame, filestem: str,
               write_fn: t.Callable[[pd.DataFrame, str], None]) -> None:
    """
    Writes the data from a DataFrame to one or more files based on the given format and specified write function.

    Parameters:
    df (pd.DataFrame): The input DataFrame containing the data to be written.
    filestem (str): A filename format string that may include placeholders for column values.
    write_fn (Callable[[pd.DataFrame, str], None]): A callable function responsible for writing each portion of the DataFrame to a file,
    taking the specific DataFrame segment and file path as arguments.

    Behavior:
    - If the `filestem` contains placeholders (keys enclosed in `{}`) that match column names in the DataFrame:
        - Validates that the DataFrame is non-empty and that all placeholder keys exist in the DataFrame columns.
        - Groups the DataFrame by the placeholder keys and writes each group to a separate file.
          The filename for each group is generated by replacing the placeholders in `filestem` with the corresponding group values.
        - The columns used for creating the groups are excluded from the output file, ensuring only the remaining columns are written.
    - If the `filestem` does not contain placeholders:
        - Writes the entire DataFrame to a single file using the specified `filestem` as the filename.
    - Creates the necessary directories for the output files if they do not already exist.
    """
    repl_keys = utils.formatkeys(filestem)
    if repl_keys:
        logging.debug(f'df columns: {df.columns}')
        logging.debug(f'df indices: {df.index.names}')
        assert len(df) != 0
        assert all([k in df.columns for k in repl_keys])

        for group, df_group in df.groupby(by=repl_keys):
            repl_vals = (group,) if isinstance(group, str) else group
            repl = dict(zip(repl_keys, repl_vals))

            filename = filestem.format(**repl)
            logging.info(f"Writing file: {filename}")
            os.makedirs(os.path.dirname(filename), exist_ok=True)

            out_cols = [c for c in df_group.columns if c not in repl_keys]

            write_fn(
                df_group[out_cols],
                filename
            )

    else:
        filename = filestem
        logging.info(f"Writing file: {filename}")
        os.makedirs(os.path.dirname(filename), exist_ok=True)

        write_fn(df, filestem)


def write_dict(df: pd.DataFrame, filestem: str, dict_depth: int) -> None:
    """
    Writes a pandas DataFrame to a dictionary file format.

    Parameters:
    df: pd.DataFrame
        The input pandas DataFrame to be written.
    filestem: str
        The base name for the output file.
    dict_depth: int
        The depth of the dictionary structure to be created.

    See frame_to_dict for details on conversion from DataFrame to dictionary.
    """
    def writeconvert(data, fname):
        np.save(fname, frame_to_dict(data, dict_depth))

    write_data(df, filestem, write_fn=writeconvert)


def write_frame(df: pd.DataFrame, filestem: str) -> None:
    """
    Writes a DataFrame to a file with a specified filestem using a custom write function.

    Parameters:
    df : pd.DataFrame
        The DataFrame to be written to a file.
    filestem : str
        The base name for the output file.

    See write_data for details.
    """
    write_data(df, filestem,
               write_fn=lambda data, fname: data.to_hdf(fname, key='corr', mode='w'))


def main(**kwargs):
    """
    Main function to handle configuration and data loading.

    Parameters:
    kwargs: Arbitrary keyword arguments
      - logging_level: str (optional)
        Optional logging level to set for the logger.
      - load_files: parameter key used to fetch configuration data.

    Behavior:
    If `kwargs` are provided:
      - Extracts logging level and attempts to fetch the DataIO configuration using the `load_files` key.

    If `kwargs` are not provided:
      - Reads the `params.yaml` file to obtain `load_files` key data to build DataIO configuration.
      - Raises a ValueError if the `load_files` key is missing.

    Returns:
      Configuration data loaded using the load function.
    """
    logging_level: str
    if kwargs:
        logging_level = kwargs.pop('logging_level', 'INFO')
        dataio_config = config.get_dataio_config(kwargs['load_files'])
    else:
        try:
            params = utils.load_param('params.yaml')['load_files']
        except KeyError:
            raise ValueError("Expecting `load_files` key in params.yaml file.")

        logging_level = params.pop('logging_level', 'INFO')
        dataio_config = config.get_dataio_config(params)

    logging.getLogger().setLevel(logging_level)

    try:
        # Try to get the running loop (will raise RuntimeError if not in a running loop)
        loop = asyncio.get_running_loop()
    except RuntimeError:
        # If no running loop, use asyncio.run() as normal
        df = asyncio.run(load(dataio_config))

    async def load_wrapper(io_config: config.DataioConfig) -> pd.DataFrame:
        df = await load(io_config)
        actions: t.Dict = dataio_config.actions
        df = [
            processor.execute(elem, actions=actions)
            for elem in df
        ]

        return pd.concat(df)

    if loop.is_running():
        # If inside a running loop (e.g., Jupyter), run the coroutine with `loop.create_task`
        return loop.create_task(load_wrapper(dataio_config))
    else:
        # If not running, use `loop.run_until_complete`
        df = loop.run_until_complete(load(dataio_config))

    actions: t.Dict = dataio_config.actions
    df = [
        processor.execute(elem, actions=actions)
        for elem in df
    ]

    df = pd.concat(df)

    return df

if __name__ == '__main__':
    print("Assuming python interpreter is being run in interactive mode.")
    print(("Result will be stored in `result` variable"
           " upon load file completion."))
    result = main()
    logging.info("Result of file load is now stored in `result` variable.")

# ------ End Output functions ------ #
