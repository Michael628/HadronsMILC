#! /usr/bin/env python3
import sys
import os
import yaml
import pickle
from python_scripts.nanny import todo_utils
from python_scripts.processing import utils

import numpy as np
from string import Formatter
import re
import logging
from dataclasses import dataclass, field
import copy


@dataclass
class RawDataProcessor:
    """Combines data in separated files into single
    nested dictionary pickle file

    Attributes
    ----------
    input : dict
        Holds `input['filestem']` and, optionally, `input['datapaths']
        used to load data from files. `filstem` can hold keyword replacements
        and `datapaths` is required for hdf5 files, providing the datapath
        to be traversed when the file is loaded
    outfilestem : str
        The output file version of `input['filestem']`
    replacements : dict
        keyword replacements that will be traversed when searching for input
        files and writing to output files
    regex : dict, optional
        regex expressions used as keyword replacements in `input['filestem']`
        to match when searching through file system
    overwrite : bool, optional
        When false and an existing output file is found, skips input files
        that correspond to data already in existing output file
    """
    input: dict
    outfilestem: dict
    replacements: dict
    regex: dict = field(default_factory=dict)
    overwrite: bool = False

    def __post_init__(self):

        # Process `replacements` down to keys needed by
        # input and output filestems
        repl_temp: dict = {}
        repl_keyset: set = {
            k
            for k in utils.formatkeys(self.input['filestem'])
        }
        repl_keyset.union({k for k in utils.formatkeys(self.outfilestem)})
        repl_temp: dict = {
            k: copy.deepcopy(v)
            for k, v in self.replacements.items()
            if k in repl_keyset
        }
        self.replacements = repl_temp

        def keysort(x):
            try:
                return ['series', 'cfg', 'time'].index(x)
            except ValueError:
                return -1

        self.input_parser = utils.FilestemFormatParser(
            filestem=self.input['filestem'],
            params=self.replacements,
            regex=self.regex,
            keysort=keysort
        )

    @property
    def writefilestem(self):
        return f"{self.outfilestem}_raw.p"

    def loadfile(self, replacements: dict = field(default_factory=dict)):

        corr = {}

        writefile: str = self.writefilestem.format(
            **dict(replacements)
        )

        if os.path.exists(writefile):

            logging.info(
                f"Loading existing output file: {writefile}"
            )

            with open(writefile, 'rb') as f:
                corr = pickle.load(f)
        else:
            logging.warning(
                f"No existing output file found at: {writefile}"
            )
            corr = {}

        return corr

    # Overwritten by subclass
    def postprocess(self, corr):
        pass

    def writefile(self, file_reps: dict,
                  corr: np.ndarray, ext: str = '') -> None:

        outfile = self.writefilestem.format(**file_reps)

        if not os.path.exists(os.path.dirname(outfile)):
            os.makedirs(os.path.dirname(outfile))

        logging.info(f"Saving file: {outfile}")

        with open(outfile, 'wb') as f:
            pickle.dump(corr, f)

    def readdata(self, corr: dict, file: str, corr_repl: dict,
                 datapaths: list[str] = None, overwrite: bool = False) -> None:
        utils.setdictval(
            corr,
            list(corr_repl.values()),
            value=utils.extractdata(file, datapaths),
            overwrite=overwrite)

    def process(self):

        for (outfile_repl, _path) in self.input_parser.traverse_replacements():

            corr = self.loadfile(outfile_repl)

            for regex_repl, infile in self.input_parser.traverse_regex():

                logging.debug(f"Processing file: {infile}")

                datapaths = self.input.get('datapaths', None)

                self.readdata(
                    corr,
                    infile,
                    regex_repl,
                    datapaths,
                    self.overwrite
                )

            if len(corr) != 0:
                self.postprocess(corr)
                self.writefile(outfile_repl, corr)


class DefaultDataProcessor(RawDataProcessor):

    @property
    def writefilestem(self):
        return self.outfilestem + "_dict.p"

    def postprocess(self, corr):

        for key, val in corr.items():
            if isinstance(corr[key], dict):
                corr[key] = utils.dict_to_corr(corr[key])

    def readdata(self, corr: dict, file: str, corr_repl: dict,
                 datapaths: list[str] = None, overwrite: bool = False) -> None:

        series_cfg = "{series}.{cfg}".format(**corr_repl)

        if 'time' not in corr_repl:
            if series_cfg not in corr or self.overwrite:
                corr[series_cfg] = utils.extractdata(file, datapaths)
        else:
            time = int(corr_repl['time'])
            if series_cfg in corr and isinstance(corr[series_cfg], dict):
                corr[series_cfg][time] = utils.extractdata(file, datapaths)
            elif series_cfg not in corr or self.overwrite:
                corr[series_cfg] = {
                    time: utils.extractdata(file, datapaths)
                }


def main():
    params = todo_utils.load_param('params.yaml')

    for run_key in params['processing']['run']:

        input_param = params['processing'][run_key]['input']
        output_param = params['processing'][run_key]['output']

        proc_param = {}
        try:
            for k, v in params['processing']['default'].items():
                proc_param[k] = copy.deepcopy(v)
        except KeyError:
            proc_param = {}

        # Grab input/output params from other run keys if requested
        if isinstance(input_param, str):
            input_param = params['processing'][input_param]['input']

        if output_param in params['processing']:
            output_param = params['processing'][output_param]['output']

        proc_param.update(params['processing'][run_key])

        data_proc: RawDataProcessor = DefaultDataProcessor(
            input=input_param,
            outfilestem=output_param,
            replacements=proc_param,
            regex=proc_param['regex'],
            overwrite=params.get('overwrite', False)
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

        data_proc.process()


if __name__ == "__main__":
    main()
