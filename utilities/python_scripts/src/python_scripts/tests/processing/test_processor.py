import pytest
import unittest.mock as mock
import processing.processor as proc
import nanny.todo_utils as todo_utils
from pprint import pprint
import pickle

@pytest.fixture
def setup_params():
    return todo_utils.load_param("../test_params/processor_tests.yaml")


def test_main(setup_params):
    print("test_main")

    setup_params['processing']['run'] = ['test']

    with mock.patch("processing.processor.todo_utils.load_param") as mock_params:
        mock_params.return_value = setup_params

        proc.main()

        with open("".join((
            "e2000n1dt2/correlators/python/m000569/",
            "vec_local/ama/corr_vec_local_ama_m000569_dict.p"
        )), 'rb') as f:
            a = pickle.load(f)
        while True:
            print([k for k in a.keys()])
            a = next(iter(a.values()))
            if not isinstance(a, dict):
                break


if __name__ == "__main__":
    pytest.main()
