import python_scripts.nanny.xml_templates.generate_lmi_params as params
import shutil
import pytest
import unittest.mock as mock
from python_scripts.nanny.xml_templates.generate_lmi_params import build_params
from python_scripts.nanny.todo_utils import load_param


def test_main():
    print("test_main")

    params = load_param('../data/test_params/hadrons_xml_test.yaml')

    return build_params(**params['lmi_param'],SERIES='a',CFG='1146')


if __name__ == "__main__":
    pytest.main()
