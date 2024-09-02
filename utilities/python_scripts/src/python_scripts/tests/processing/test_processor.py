import pytest
import unittest.mock as mock
import processing.processor as proc
import nanny.todo_utils as todo_utils
from pprint import pprint


@pytest.fixture
def setup_params():
    return todo_utils.load_param("../../processing/example_params.yaml")


def test_main(setup_params):
    print("test_main")

    params = setup_params

    with mock.patch("processing.processor.todo_utils.load_param") as mock_params:
        mock_params.return_value = params

        proc.main()


if __name__ == "__main__":
    pytest.main()
