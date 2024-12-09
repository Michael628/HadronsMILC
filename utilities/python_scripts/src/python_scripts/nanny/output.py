from dataclasses import dataclass
from python_scripts.config import ConfigBase
import typing as t
from python_scripts import utils


@dataclass
class FileOutputConfig(ConfigBase):


def create_epack_output(files: t.Dict):
    pass

def create_meson_output():
    pass

def create_high_modes_output():
    pass


def get_output(output_label: str):
    output = {
        "epack": create_epack_output,
        "meson": create_meson_output,
        "high_modes": create_high_modes_output,
    }

    if output_label in output:
        return output[output_label]
    else:
        raise ValueError(f"No config implementation for `{output_label}`.")


def get_epack_outfiles(tasks: config.ConfigBase, file_params: t.Dict,
                       run_config: ConfigBase):

    home = file_params['home']

    eigfile_params: t.Dict
    if tasks.multifile:
        eigfile_params = file_params['epack']['eig']
    else:
        eigfile_params = file_params['epack']['eigdir']

    replacements: t.Dict
    goodsize: int
    filestem: str

    replacements = eigfile_params.get('replacements', {})
    goodsize = eigfile_params['goodize']
    filestem = eigfile_params['filestem']

    replacements = utils.process_params(**replacements)
    replacements = run_config.__dict__.update(replacements)

    utils.process_files()
    epack_tasks = tasks['epack']
    file_params = param['files']['epack']

    generate_eigs: bool = not epack_tasks['load']

    multifile: bool
    multifile = epack_tasks.get('multifile', False)

    save_evals: bool
    save_evals = generate_eigs or epack_tasks.get('save_evals', False)

    files = []
    if generate_eigs:
        eig_key = 'eigdir' if multifile else 'eig'
        files.append(file_params[eig_key])
        files.append(file_params['eval'])

        output = dict(files[key].items())
        if multifile:
            output['formatting'] = {
                'eig_index': [str(i) for i in range()]
            }
    if save_evals:
        ret.append(files['eval'])

    return ret

def get_outfile_config(file_label: str):

    output = {
        "epack": create_epack_output,
        "meson": create_meson_output,
        "high_modes": create_high_modes_output,
    }

    if output_label in output:
        return output[output_label]
    else:
        raise ValueError(f"No config implementation for `{output_label}`.")

