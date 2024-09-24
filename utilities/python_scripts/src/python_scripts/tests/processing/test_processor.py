import shutil
import pytest
import unittest.mock as mock
import processing.processor as proc
import nanny.todo_utils as todo_utils
import python_scripts

from pprint import pprint
import pickle
import numpy as np


@pytest.fixture
def setup_params():
    return todo_utils.load_param("../data/test_params/processor_tests.yaml")


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
    pytest.main()
