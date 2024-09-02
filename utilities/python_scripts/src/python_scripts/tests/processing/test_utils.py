import pytest
import unittest.mock as mock
import processing.utils as utils
import nanny.todo_utils as todo_utils
from pprint import pprint


@pytest.fixture
def setup_params():
    return todo_utils.load_param("../../processing/example_params.yaml")


def test_dictval_iter(setup_params):

    print("test_dictval_iter")

    params = setup_params['processing']['test']

    params.update(params['regex'])

    def keysort(val):
        if 'key_order' in params and val in params['key_order']:
            return params['key_order'].index(val)
        else:
            return -1

    sub_dict = utils.formatdict(params['input']['filestem'], params)

    keys, val_iter = utils.dictval_iter(sub_dict,key_sort=keysort)

    for vals in val_iter:
        for k, v in dict(zip(keys, vals)).items():
            print(k)




if __name__ == "__main__":
    pytest.main()
