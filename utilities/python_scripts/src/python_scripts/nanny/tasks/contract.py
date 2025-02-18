import logging
import os
import typing as t
from dataclasses import dataclass

import pandas as pd

from python_scripts import config as c, utils
from python_scripts.nanny import SubmitConfig, TaskBase
from python_scripts.nanny.config import OutfileList
from python_scripts.a2a.config import DiagramConfig


@c.dataclass_with_getters
class SubmitContractConfig(SubmitConfig):
    _diagram_params: t.Dict[str, DiagramConfig]
    hardware: t.Optional[str] = None
    logging_level: t.Optional[str] = None

    def __init__(self, **kwargs):
        self.hardware = kwargs.get('hardware', None)
        self.logging_level = kwargs.get('logging_level', None)
        self._diagram_params = {}
        for k, v in kwargs['_diagram_params'].items():
            self._diagram_params[k] = DiagramConfig.create(**v)


@dataclass
class ContractTask(TaskBase):
    diagrams: t.List[str]


def input_params(task_config: ContractTask, submit_config: SubmitContractConfig,
                 outfile_config_list: OutfileList) -> t.Tuple[t.List[str], None]:
    input_yaml = submit_config.public_dict
    input_yaml['diagrams'] = {}
    for diagram in task_config.diagrams:
        d_params = submit_config.diagram_params[diagram]
        d_params.set_filenames(outfile_config_list)
        input_yaml['diagrams'][diagram] = d_params.string_dict()

    return input_yaml, None


def catalog_files(task_config: ContractTask, submit_config: SubmitContractConfig,
                  outfile_config_list: OutfileList) -> t.List[str]:
    def build_row(filepath: str, repls: t.Dict[str, str]) -> t.Dict[str, str]:
        repls['filepath'] = filepath
        return repls

    replacements = submit_config.string_dict()
    outfile_config = outfile_config_list.contract
    df = []
    outfile = outfile_config.filestem + outfile_config.ext
    filekeys = utils.formatkeys(outfile)
    for diagram in task_config.diagrams:
        d_params = submit_config.diagram_params[diagram]
        d_params.set_filenames(outfile_config_list)
        replacements['mass'] = d_params.mass
        replacements['gamma'] = d_params.diagram_label
        files = utils.process_files(
            outfile,
            processor=build_row,
            replacements={k: v for k, v in replacements.items() if k in filekeys}
        )
        dict_of_rows = {k: [file[k] for file in files] for k in files[0] if len(files) > 0}

        new_df = pd.DataFrame(dict_of_rows)
        new_df['good_size'] = outfile_config.good_size
        new_df['exists'] = new_df['filepath'].apply(os.path.exists)
        new_df['file_size'] = None
        new_df.loc[new_df['exists'], 'file_size'] = new_df[new_df['exists']]['filepath'].apply(os.path.getsize)
        df.append(new_df)

    df = pd.concat(df, ignore_index=True)

    return df


def bad_files(task_config: ContractTask, submit_config: SubmitContractConfig,
              outfile_config_list: OutfileList) -> t.List[str]:
    df = catalog_files(submit_config, task_config, outfile_config_list)

    return list(df[(df['file_size'] >= df['good_size']) != True]['filepath'])


def get_task_factory():
    return ContractTask


def get_submit_factory() -> t.Callable[..., SubmitContractConfig]:
    return SubmitContractConfig.create


def processing_params(task_config: ContractTask, submit_config: SubmitContractConfig,
                      outfile_config_list: OutfileList) -> t.Dict:
    infile_stem = outfile_config_list.contract.filename
    infile_stem = infile_stem.replace("{gamma}","{diagram_label}")
    outfile = outfile_config_list.contract.filestem
    outfile = outfile.replace("{gamma}","{diagram_label}")
    filekeys = utils.formatkeys(infile_stem)
    proc_params = {'run': task_config.diagrams}
    outfile = outfile.replace("correlators", "dataframes")
    outfile = outfile.replace("_{series}", "")
    outfile += ".h5"
    replacements = {k: v for k, v in submit_config.string_dict().items() if k in filekeys}

    for k in task_config.diagrams:
        diagram_dict = submit_config.diagram_params[k].string_dict()
        replacements.update({k: v for k, v in diagram_dict.items() if k in filekeys})
        proc_params[k] = {
            "logging_level": getattr(submit_config, "logging_level", "INFO"),
            "load_files": {
                "filestem": infile_stem,
                "regex": {
                    "series": "[a-z]",
                    "cfg": "[0-9]+"
                },
                "dict_labels": ["seedkey", "gamma"],
            },
            "actions": {
                "drop": "seedkey"
            },
            "out_files": {
                "filestem": outfile,
                "type": "dataframe"
            },
        }

        proc_params[k]['load_files']["replacements"] = replacements.copy()
        proc_params[k]['load_files']["replacements"]["diagram_label"] = diagram_dict['diagram_label']
        proc_params[k]['load_files']["array_params"] = {
            "order": ["t"],
            "labels": {
                "t": "0..47"
            }
        }

    return proc_params
