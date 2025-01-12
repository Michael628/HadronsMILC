from python_scripts.processing.dataio import DataLoad
import python_scripts
import python_scripts.utils as utils
import copy
import logging

import typing as t


def main(param_input: t.Union[t.Dict,str]):
    try:
        params = utils.load_param(param_input)
    except TypeError:
        params = param_input

    res = {}
    for run_key in params['processing']['run']:

        run_params = params['processing'][run_key]

        input_params = run_params.pop('input')

        data_keys = input_params.pop('data_keys')
        data_labels = input_params.pop('data_labels')
        datapaths = input_params.pop('datapaths')

        # try:
        # input_param = params['processing'][input_param]['input']
        # except TypeError:
        # pass

        if 'default' in params['processing']:
            input_params.update((
                (k, copy.deepcopy(v))
                for k, v in params['processing']['default'].items()
                if k not in input_params
            ))

        loader: DataLoad = DataLoad(
            input_params=input_params,
            data_keys=data_keys,
            data_labels=data_labels,
            datapaths=datapaths
        )

        python_scripts.setup()

        logging_level = params['processing'].get('logging_level', logging.INFO)


        res[run_key] = loader.execute(,

    return res


if __name__ == "__main__":
    main('params.yaml')
