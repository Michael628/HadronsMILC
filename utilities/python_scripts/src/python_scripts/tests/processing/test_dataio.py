import logging
import os
import sys
import python_scripts.processing.dataio as dataio
import pdb
from pprint import pprint
import pickle
import numpy as np
import pandas as pd
from python_scripts.processing.dataio import (
    string_replacement_gen,
    file_regex_gen,
    ndarray_to_frame,
    load_pickle,
    load_input
)

logging.basicConfig(
    format="%(asctime)s - %(levelname)-5s - %(message)s",
    style="%",
    datefmt="%Y-%m-%d %H:%M:%S",
    level="DEBUG",
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)


def create_directory_for_file(file_path):
    # Get the directory name from the file path
    directory = os.path.dirname(file_path)

    # Create the directory if it does not exist
    if not os.path.exists(directory):
        os.makedirs(directory)  # Creates all intermediate directories


def test_replacements():
    print("test_replacements")

    params = {
        'input_params': {
            'filestem': ("/home/michael/perlmutter_sync/qed/a2a_corrs/sib/"
                         "corr_sib_{block}_{series}.{cfg}.p"),
            'replacements': {
                # ,'LLLHHL','LHHLLL','HLLLLH']
                'block': ['LHHLLL_seed20', 'LLLLLL'],
                # 'series': 'a',
                # 'cfg': '60'
            },
            'regex': {
                'series': '[a-z]',
                'cfg': '[0-9]+'
                # 'cfg': '122[0-9]'
            },
        },
        'preprocess': lambda x: x[0, 0],
        'data_keys': [
            'seedkey',
            'gamma',
            'dt'
        ],
        'data_labels': {
            'dt': '0..47'
        }
    }

    for a, b in string_replacement_gen(
            params['input_params']['filestem'],
            params['input_params']['replacements']):

        for c, d in file_regex_gen(b, params['input_params']['regex']):
            pprint(d)
            pprint(a)
            pprint(c)
    # pdb.set_trace()
    # dl = dataio.DataLoad(**params)
    # df: pd.DataFrame = dl.execute()

    # df[['re','im']] = df['corr'].apply(lambda x: pd.Series([x.real,x.imag]))
    # df.drop(columns=['corr'],inplace=True)
    # pprint(df)

    # outfile = "../data/sib/corr_sib.feather"
    # create_directory_for_file(outfile)
    # df.to_feather(outfile)


def test_pickle_to_df():
    print("test_pickle_to_df")

    params = {
        'input_params': {
            'filestem': ("/home/michael/perlmutter_sync/qed/a2a_corrs/sib/"
                         "corr_sib_LLLLLL_a.60.p")
        },
        'preprocess': lambda x: x[0, 0],
        'data_keys': [
            'seedkey',
            'gamma',
            'dt'
        ],
        'data_labels': {
            'dt': '0..47'
        }
    }

    with open(params['input_params']['filestem'], 'rb') as f:
        a = pickle.load(f)

    b = ndarray_to_frame(
        next(iter(a.values()))['GX_GX'],
        label_order=['t1', 't2', 'dt'],
        array_labels={
            't1': '0..47',
            't2': '0..47',
            'dt': '0..47'
        })

    pprint(b._data)


def test_load_pickle():
    print("test_pickle_to_df")

    params = {
        'input_params': {
            'filestem': ("/home/michael/perlmutter_sync/qed/a2a_corrs/sib/"
                         "corr_sib_LLLLLL_a.60.p")
        },
        # 'preprocess': lambda x: x[0, 0],
        'data_order': [
            'seedkey',
            'gamma',
            't1',
            't2',
            'dt'
        ],
        'array_labels': {
            't1': '0..47',
            't2': '0..47',
            'dt': '0..47'
        }
    }

    return load_pickle(
        filename=params['input_params']['filestem'],
        data_order=params['data_order'],
        # preprocess=params['preprocess'],
        array_labels=params['array_labels']
    )


def test_load_input():
    print("test_pickle_to_df")

    params = {
        'input_params': {
            'filestem': ("/home/michael/perlmutter_sync/qed/a2a_corrs/sib/"
                         "corr_sib_LLLLLL_a.60.p")
        },
        # 'preprocess': lambda x: x[0, 0],
        'data_order': [
            'seedkey',
            'gamma',
            't1',
            't2',
            'dt'
        ],
        'array_labels': {
            't1': '0..47',
            't2': '0..47',
            'dt': '0..47'
        }
    }

    df = load_input(config=params.pop('input_params'), **params)

    return df


if __name__ == "__main__":
    test_extract_pickle()
