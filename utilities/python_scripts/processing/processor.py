#! /usr/bin/env python3
import sys
import os
import yaml
import pickle
import nanny.todo_utils as todo_utils
import processing.utils as utils
import numpy as np
from string import Formatter
import re
import logging
from dataclasses import dataclass, field
import copy


@dataclass
class FilestemFormatParser:
    filestem: str
    params: dict = field(default_factory=dict)
    keysort: callable = None

    def __post_init__(self):
        self.filekeys = utils.formatkeys(self.filestem, key_sort=self.keysort)

    @property
    def keys(self):
        return self.filekeys

    def filestem_keyorder(self, filename_only: bool = False) -> dict:
        """Returns a list of keys in the order they are found in `filestem`

        Parameters
        ----------
        filename_only : bool, optional
            If true, `filestem` is first split using `os.path.basename`
            to get only the filename.
        """
        if filename_only:
            filestr = os.path.basename(self.filestem)
        else:
            filestr = self.filestem

        order = [
            k
            for k in self.filekeys
            if k in filestr
        ]
        order.sort(key=filestr.index)

        return order

    def format(self, repl: dict) -> str:
        return self.filestem.format(**repl)

    def traverse_replacements(self) -> (str, dict):
        _, val_iter = utils.format_iter(
            self.filestem,
            self.params,
            self.keysort
        )

        for vals in val_iter:
            repl = dict(zip(self.filekeys, vals))
            yield self.format(repl), repl


@dataclass
class RawDataProcessor:
    input: dict
    outfiles: dict
    params: dict
    regex: dict = field(default_factory=dict)
    overwrite: bool = False
    keysort: callable = None

    def __post_init__(self):

        self.check_output_keys()

        # self.output_parsers = {
        #     key: FilestemFormatParser(val, self.params)
        #     for key, val in self.output.items()
        # }

        input_params = copy.deepcopy(self.params)

        # Add parentheses to regex elements if not already present
        for key, val in self.regex.items():
            if val[0] != '(' or val[-1] != ')':
                self.regex[key] = f"({val})"

        input_params.update(self.regex)

        self.input_parser = FilestemFormatParser(
            self.input['filestem'],
            input_params,
            self.keysort
        )

        output_keys = utils.formatkeys(next(iter(self.outfiles.values())))
        input_keys = self.input_parser.keys

        self.data_keys = [k for k in input_keys if k not in output_keys]

        if any(k not in self.regex for k in self.data_keys):
            raise KeyError(
                "Expecting regex replacements for keys: {}".format(
                    ", ".join(self.data_keys)
                )
            )

    def check_output_keys(self) -> None:
        """Ensure that all format strings in `output` have
        the same replacement variables
        """

        key_set = {
            " ".join(utils.formatkeys(val))
            for val in self.outfiles.values()
        }

        if len(key_set) != 1:
            raise ValueError(" ".join((
                "All output string replacements",
                "must be the same.")))

    # def postprocess(self, corr: np.ndarray) -> None:
    #     if 'time' in index_dict:
    #         time = int(re_match[index_dict['time']])
    #         temp = utils.extractdata(
    #             f"{directory}/{file}",
    #             datapaths
    #         )

    #         if seriescfg not in corr:
    #             corr[seriescfg] = np.zeros((temp.shape[-1],)+temp.shape,dtype=np.complex128)

    #         corr[seriescfg][time] = temp
    #     else:
    #         corr[seriescfg] = utils.extractdata(f"{directory.format(**file_reps)}/{file}",datapaths)

    def writefile(self, file_reps: dict, corr: np.ndarray) -> None:

        # self.postprocess(corr)

        for outtype, outfilestem in self.outfiles.items():
            outfile = outfilestem.format(**file_reps)

            print(outtype,outfile)

            if not os.path.exists(os.path.dirname(outfile)):
                os.makedirs(os.path.dirname(outfile))

            logging.info(f"Writing file: {outfile}")
            with open(outfile,'wb') as f:
                if outtype == 'numpy':
                    pickle.dump(utils.dict_to_corr(corr), f)
                elif outtype == 'dict':
                    pickle.dump(corr, f)

    def readdata(self, corr: np.ndarray, file: str, regex_repl: dict,
                 datapaths: list[str] = None, overwrite: bool = False) -> None:
        utils.setdictval(
            corr,
            list(regex_repl.values()),
            value=utils.extractdata(file, datapaths),
            overwrite=overwrite)

    def process(self):

        replacements: dict = {}

        for infile_path, infile_repl in self.input_parser.traverse_replacements():

            print(f"infile_path:{infile_path}")
            file_complete: bool = len(replacements) == 0 or any(
                replacements[k] != infile_repl[k]
                for k in replacements.keys()
            )

            if file_complete:
                print("file complete")
                if len(replacements) != 0:
                    print("writing file!")
                    self.writefile(replacements, corr)

                print(f"replacements: {replacements}")
                replacements: dict = infile_repl
                corr: dict = {}

                # Check for existing output file
                if 'dict' in self.outfiles:
                    dict_outfile: str = self.outfiles['dict'].format(
                        **infile_repl
                        )

                    if os.path.exists(dict_outfile) and not self.overwrite:

                        logging.info(
                            f"Loading existing dictionary: {dict_outfile}"
                        )

                        with open(dict_outfile, 'rb') as f:
                            corr = pickle.load(f)

            # FIX ME: Assumes all regex matches occur in file name,
            # not in the directory path.
            infile_directory, infile_match = os.path.split(infile_path)

            pattern = re.compile(infile_match)

            files: list[str] = os.listdir(infile_directory)

            regex_keyorder = self.input_parser.filestem_keyorder(True)

            _ = [
                regex_keyorder.remove(k)
                for k in regex_keyorder
                if k not in self.data_keys
            ]

            for file in files:
                re_match = pattern.match(file)

                if re_match is not None:

                    logging.debug(f"Processing file: {file}")

                    regex_repl = dict(
                        (k, re_match[regex_keyorder.index(k)])
                        for k in self.data_keys)

                    datapaths = self.input.get('datapaths', None)

                    self.readdata(corr, f"{infile_directory}/{file}", regex_repl,
                                  datapaths, self.overwrite)


def main():
    params = todo_utils.load_param('params.yaml')

    for run_key in params['processing']['run']:

        input_param = params['processing'][run_key]['input']
        output_param = params['processing'][run_key]['output']

        try:
            proc_param = params['processing']['default']
        except KeyError:
            proc_param = {}

        # Grab input/output params from other run keys if requested
        if isinstance(input_param, str):
            input_param = params['processing'][input_param]['input']
        if isinstance(output_param, str):
            output_param = params['processing'][output_param]['output']

        proc_param.update(params['processing'][run_key])

        def keysort(val):
            if 'key_order' in proc_param and val in proc_param['key_order']:
                return proc_param['key_order'].index(val)
            else:
                return -1

        rawdata_proc = RawDataProcessor(
            input=input_param,
            outfiles=output_param,
            regex=proc_param['regex'],
            overwrite='overwrite' in params and params['overwrite'],
            keysort=keysort,
            params=proc_param
        )

        if 'logging_level' in params['processing']:
            logging_level = params['processing']['logging_level']
        else:
            logging_level = logging.INFO

        logging.basicConfig(
            format="%(asctime)s - %(levelname)-5s - %(message)s",
            style="%",
            datefmt="%Y-%m-%d %H:%M:%S",
            level=logging_level,
            handlers=[
                logging.StreamHandler(sys.stdout)
            ]
        )

        rawdata_proc.process()


if __name__ == "__main__":
    main()
