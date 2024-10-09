import python_scripts.processing.dataio as dataio
import pdb
from pprint import pprint
import pickle
import numpy as np


def test_main():
    print("test_main")

    params = {
        'input_params': {
            'filestem': ("/home/michael/perlmutter_sync/qed/a2a_corrs/sib/"
                         "corr_sib_LLLLLL_a.60.p"),
            # 'replacements': {
            # 'block':['LHHLLL','HLLHHHH','HHHLLH'],#,'LLLHHL','LHHLLL','HLLLLH']
            # 'series': 'a',
            # 'cfg': '60'
            # },
            # 'regex': {
                # 'series': '[a-z]',
                # 'cfg': '[0-9]+'
                # 'cfg': '122[0-9]'
            # },
        },
        # 'preprocess': time_average,
        'data_keys': [
            'seedkey',
            'gamma',
        ],
    }
    pdb.set_trace()
    dl = dataio.DataLoad(**params)
    df = dl.execute()
    pprint(df)


if __name__ == "__main__":
    test_main()
