# This file generates xml parameters for the HadronsMILC app.
# Tasks performed:
#
# 1: Load eigenvectors
# 2: Generate noise sources
# 3: Solve low-mode propagation of sources
# 4: Solve CG on result of step 3
# 5: Subtract 3 from 4
# 6: Save result of 5 to disk

import logging
import typing as t

from pydantic.dataclasses import dataclass

from python_scripts.nanny import TaskBase
from python_scripts.nanny.config import OutfileList
from python_scripts.nanny.tasks.hadrons import SubmitHadronsConfig, templates


@dataclass
class TestTask(TaskBase):
    mass: str
    gammas: str


def input_params(
    tasks: TestTask,
    submit_config: SubmitHadronsConfig,
    outfile_config_list: OutfileList,
) -> t.Tuple[t.List[t.Dict], t.Optional[t.List[str]]]:
    submit_conf_dict = submit_config.string_dict()

    modules = []

    meson_path = outfile_config_list.meson_ll.filestem

    nvecs = str(3 * submit_config.time)

    for seed_index in range(submit_config.noise):
        w_name = f"w{seed_index}"
        modules.append(templates.time_diluted_noise(w_name, 1))
        mass_label = tasks.mass

        outfile = meson_path.format(
            mass=submit_config.mass_out_label[mass_label],
            w_index="test",
            v_index="test",
            **submit_conf_dict,
        )

        modules.append(
            templates.meson_field(
                name=f"mf_{seed_index}_{seed_index}",
                action=f"stag_mass_{mass_label}",
                block=nvecs,
                gammas=tasks.gammas,
                gauge="",
                low_modes="",
                left=w_name + "_vec",
                right=w_name + "_vec",
                output=outfile,
                apply_g5="false",
            )
        )

    return modules, None


def bad_files(
    task_config: TaskBase,
    submit_config: SubmitHadronsConfig,
    outfile_config_list: OutfileList,
) -> t.List[str]:
    logging.warning(
        "Check completion succeeds automatically. No implementation of bad_files function in `hadrons_a2a_vectors.py`."
    )
    return []


def get_task_factory():
    return TestTask
