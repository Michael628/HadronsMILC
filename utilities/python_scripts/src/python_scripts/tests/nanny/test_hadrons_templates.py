import shutil
import pytest
import unittest.mock as mock
from python_scripts.nanny.xml_templates.generate_lmi import build_params
from python_scripts.nanny.todo_utils import load_param
from pprint import pprint
from dict2xml import dict2xml as dxml

def test_main():
    print("test_main")

    params = load_param('../data/test_params/hadrons_xml_test.yaml')

    tasks = params['job_setup']['L'].get('tasks', {})

    xml_dict, schedule = build_params(
        tasks=tasks,
        SERIES='a', CFG='1146',
        **params['lmi_param']
    )

    with open('../data/test_xml/hadrons_xml_test.xml', 'w') as f:
        print(dxml(xml_dict), file=f)
    with open('../data/test_schedule/hadrons_xml_test.sched', 'w') as f:
        f.write(str(len(schedule))+"\n"+"\n".join(schedule))


if __name__ == "__main__":
    test_main()
