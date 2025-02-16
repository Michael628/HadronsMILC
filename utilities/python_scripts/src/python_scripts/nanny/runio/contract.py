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
    _diagram_params: t.Dict
    hardware: t.Optional[str] = None
    logging_level: t.Optional[str] = None


@dataclass
class ContractTask(TaskBase):
    diagrams: t.List[str]


def build_params(submit_config: SubmitContractConfig,
              task_config: ContractTask, outfile_config_list: OutfileList) -> t.Tuple[t.List[str], None]:
    input_yaml = submit_config.public_dict
    input_yaml['diagrams'] = {}
    for diagram in task_config.diagrams:
        d_params = submit_config.diagram_params[diagram]
        d = DiagramConfig.create(outfile_config_list,**d_params)
        input_yaml['diagrams'][diagram] = d.string_dict()

    return input_yaml, None

def catalog_files(submit_config: SubmitContractConfig,
              task_config: ContractTask, outfile_config_list: OutfileList) -> t.List[str]:

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
        d = DiagramConfig.create(**d_params)
        replacements['mass'] = d.mass
        replacements['gamma'] = d.diagram_label
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


def bad_files(submit_config: SubmitContractConfig,
              task_config: ContractTask, outfile_config_list: OutfileList) -> t.List[str]:
    df = catalog_files(submit_config, task_config, outfile_config_list)

    return list(df[(df['file_size'] >= df['good_size']) != True]['filepath'])


def get_task_factory():
    return ContractTask

def get_submit_factory() -> t.Callable[..., SubmitContractConfig]:
    return SubmitContractConfig.create