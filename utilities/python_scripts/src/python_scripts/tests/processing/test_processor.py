import shutil
import pytest
import unittest.mock as mock
import python_scripts.processing.processor as proc
import python_scripts.nanny.todo_utils as todo_utils
import python_scripts
import python_scripts.processing.dataio as dataio

from pprint import pprint
import pickle
import numpy as np


@pytest.fixture
def setup_params():
    return todo_utils.load_param("../data/test_params/processor_tests.yaml")


params = {
    'vec_local_ranLL_lv': {
        'input_params': {
            'filestem': ("/home/michael/perlmutter_sync/l12896/"
                         "e8000n1dt1_python/m{mass}/{diagram}/{dset}_e{eigs}/"
                         "corr_{diagram}_{dset}_m{mass}_dict.p"),
            'replacements': {
                'dset': 'ranLL',
                'mass': '001326',
                'diagram': 'vec_local',
                'eigs': ['2000', '4000', '8000'],
            },
        },
        'data_order': [
            'series.cfg',
            'tsource',
            'gamma',
            'dt'
        ],
        'array_labels': {
            'tsource': '0..95',
            'gamma': ['GX_GX', 'GY_GY', 'GZ_GZ'],
            'dt': '0..95'
        }
    },
    'vec_local_a2aLL_lv': {
        'input_params': {
            'filestem': ("/home/michael/perlmutter_sync/l12896/"
                         "e8000n1dt1_python/m{mass}/{diagram}/{dset}_e{eigs}/"
                         "corr_{diagram}_{dset}_m{mass}_dict.p"),
            'replacements': {
                'dset': 'a2aLL',
                'mass': '001326',
                'diagram': 'vec_local',
                'eigs': ['2000', '4000', '8000'],
            },
        },
        'data_order': [
            'series.cfg',
            'gamma',
            'data'
        ],
        'array_labels': {
            'gamma': ['GX_GX', 'GY_GY', 'GZ_GZ']
        }
    },
    'vec_local_ama_lv': {
        'input_params': {
            'filestem': ("/home/michael/perlmutter_sync/l12896/"
                         "e{eigs}n1dt1_python/m{mass}/{diagram}/{dset}/"
                         "corr_{diagram}_{dset}_m{mass}_dict.p"),
            'replacements': {
                'dset': 'ama',
                'mass': '001326',
                'diagram': 'vec_local',
                'eigs': '8000'
            },
        },
        'data_order': [
            'series.cfg',
            'time',
            'gamma',
            'data'
        ],
        'array_labels': {
            'time': '0..95',
            'gamma': ['GX_GX', 'GY_GY', 'GZ_GZ']
        }
    },
    'vec_local_high_sv': {
        'input_params': {
            'filestem': ("/home/michael/lattice_data/callat/correlators/"
                         "{diagram}/{dset}/"
                         "{dset}_{mass}_{diagram}_{cfgset}.p"),
            'replacements': {
                'diagram': 'local_vec',
                'cfgset': [
                    "b.6_3090_dconf24", "b.18_3078_dconf24",
                    "a.6_3102_dconf24", "a.18_3090_dconf24",
                    "a.3108_3516_dconf12", "a.3528_6108_dconf12",
                    "b.3108_6108_dconf12"
                ],
                'mass': "m001326",
                'dset': ['ama', 'ranLL'],
            }
        },
        'data_order': [
            'series.cfg',
            'gamma',
            'time',
            'data'
        ],
    },
    'vec_local_a2aLL_sv': {
        'input_params': {
            'filestem': ("/home/michael/lattice_data/callat/correlators/"
                         "{diagram}/{dset}/"
                         "{dset}_{mass}_{diagram}_{cfgset}.p"),
            'replacements': {
                'diagram': 'local_vec',
                'cfgset': [
                    "b.6_3090_dconf24", "b.18_3078_dconf24",
                    "a.6_3102_dconf24", "a.18_3090_dconf24",
                    "a.3108_3516_dconf12", "a.3528_6108_dconf12",
                    "b.3108_6108_dconf12"
                ],
                'mass': "m001326",
                'dset': ['a2aLL'],
            }
        },
        'data_order': [
            'series.cfg',
            'gamma',
            'data'
        ],
    },
}


def test_processor():
    input_params = params['vec_local_ranLL_lv'].pop('input_params')

    return dataio.load_input(input_params, **params['vec_local_ranLL_lv'])


def test_main(setup_params):
    print("test_main")

    setup_params['processing']['run'] = ['test']

    params = setup_params['processing']['test']

    params.update(setup_params['processing']['default'])

    def keysort(x):
        try:
            return ['series', 'cfg', 'time'].index(x)
        except ValueError:
            return -1

    # input_parser = FileParser(
    #     filestem=params['input']['filestem'],
    #     replacements=params,
    #     regex=params['regex'],
    #     keysort=keysort
    # )

    load_str = "processing.processor.todo_utils.load_param"
    with mock.patch(load_str) as mock_params:
        mock_params.return_value = setup_params

        try:
            shutil.rmtree("../data/e2000n1dt2")
        except OSError:
            pass

        main('params.yaml')

        # for reps, outstem in input_parser.traverse_replacements():
        #     print(outstem)


if __name__ == "__main__":
    # pytest.main()
    a = test_processor()
