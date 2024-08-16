#! /usr/bin/env python3

import yaml
import pickle
from todo_utils import load_param
from processingutils import dict_traversal_iterator
import numpy as np
from string import Formatter
import re

def extract_data(file):

    if file[-1] == "p":
        return pickle.load(open(file,"rb"))
    elif file[-2:] == "h5":
        raise Exception("hdf5 needs to be implemented.")
    else:
        raise Exception("Only pickle and hdf5 files are supported.")

def main():
    """
    Parameters
    ----------
    replacements : dict
        key-val replacements for formatting input/output file names.
        All lists are iterated over for replacement
    input : dict
        key : str
            Expected keys for input dict are ama, ranLL, a2aLL
        value : str
            The format string containing the file from which to
            read the data. Data is expected to be a dictionary keyed
            by 
    output : dict
        key : str
            Expected keys for input dict are ama, ranLL, a2aLL
        value : str
            The format string containing the file to write the
            desired data

    Example
    -------
    An example yaml parameter file would have the follow format
    processing:
      gamma: 'sib'
      mass: '002426'
      eigs: '1000'
      dt: '1'
      noise: '1'
      input:
        a2aLL: a2a_corrs/{gamma}/corr_{gamma}_LLLLLL_{series}.{cfg}.p
      output:
        a2aLL:
          numpy: e{eigs}n{noise}dt{dt}/correlators/python/{gamma}/a2a/corr_{gamma}_LLLLLL_numpy.p
          dict:  e{eigs}n{noise}dt{dt}/correlators/python/{gamma}/a2a/corr_{gamma}_LLLLLL.p

    """
    params = load_param('params.yaml')

    proc_params = params['processing']

    corr_out = {}

    input_dsets = list(proc_params['input'].keys())
    output_dsets = list(proc_params['output'].keys())


    for dset in input_dsets:
        format_keys, replacement_iter = build_traversal_list(build_format_dict(proc_params['input']))
    for replacement in replacement_iter:

        file_reps = dict(zip(format_keys,replacement))

            print(dset)

            corr_dict = {}

            # list of required formatting keys for input string
            input_keys = [k[1] for k in Formatter().parse(proc_params['input'][dset]) if k[1] is not None]

            directory, infile_match = os.path.split(proc_params['input'][dset])

            files = os.listdir(directory.format(**file_reps))

            configs = set()

            pattern = re.compile(infile_match.format(series="([a-z])",cfg="([0-9]*)",**file_reps))

            for file in files:
                re_match = pattern.match(file)

                if re_match != None:
                    series = re_match[0]
                    config = re_match[1]
                    seriescfg = f'{series}.{config}'

                    corr_dict[seriescfg] = extract_data(f"{directory.format(**file_reps)}/{file}")
            
            if dset in output_dsets:
                for outtype, outfile in proc_params['output'][dset].items():
                    if not os.path.exists(os.path.dirname(outfile.format(**file_reps))):
                        os.makedirs(os.path.dirname(outfile.format(**file_reps)))
                    if outtype == 'numpy':
                        pickle.dump(dict_to_corr(corr_dict),open(outfile.format(**file_reps),'wb'))
                    elif outtype == 'dict':
                        pickle.dump(corr_dict,open(outfile.format(**file_reps),'wb'))                    
            
if __name__ == "__main__":
    main()