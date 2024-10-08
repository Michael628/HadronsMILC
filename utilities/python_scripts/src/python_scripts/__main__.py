from python_scripts.processing.dataio import DataLoad


def main(param_input: dict | str):
    try:
        params = todo_utils.load_param(param_input)
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

        logging_level = params['processing'].get('logging_level', logging.INFO)

        logging.basicConfig(
            format="%(asctime)s - %(levelname)-5s - %(message)s",
            style="%",
            datefmt="%Y-%m-%d %H:%M:%S",
            level=logging_level,
            handlers=[
                logging.StreamHandler(sys.stdout)
            ]
        )

        res[run_key] = loader.execute()

    return res


if __name__ == "__main__":
    main('params.yaml')
