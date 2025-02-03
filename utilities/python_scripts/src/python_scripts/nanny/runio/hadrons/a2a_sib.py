# This file generates xml parameters for the HadronsMILC app.
# Tasks performed:
#
# 1: Load eigenvectors
# 2: Generate noise sources
# 3: Solve low-mode propagation of sources
# 4: Solve CG on result of step 3
# 5: Subtract 3 from 4
# 6: Save result of 5 to disk
import itertools
from dataclasses import dataclass
import logging
from python_scripts.nanny.config import OutfileList
from python_scripts.nanny.runio.hadrons import SubmitHadronsConfig, templates
from python_scripts.nanny import TaskBase
import typing as t


@dataclass
class A2ASIBTask(TaskBase):
    mass: str
    gammas: str
    multifile: bool = True
    low_memory_mode: bool = True
    w_indices: t.Optional[t.List[int]] = None
    v_indices: t.Optional[t.List[int]] = None


def build_params(submit_config: SubmitHadronsConfig, tasks: A2ASIBTask,
                 outfile_config_list: OutfileList) -> t.Tuple[t.List[t.Dict], t.Optional[t.List[str]]]:
    submit_conf_dict = submit_config.string_dict()

    gauge_filepath = outfile_config_list.gauge_links.filestem.format(**submit_conf_dict)
    gauge_fat_filepath = outfile_config_list.fat_links.filestem.format(**submit_conf_dict)
    gauge_long_filepath = outfile_config_list.long_links.filestem.format(**submit_conf_dict)

    modules = [
        templates.load_gauge('gauge', gauge_filepath),
        templates.load_gauge('gauge_fat', gauge_fat_filepath),
        templates.load_gauge('gauge_long', gauge_long_filepath)
    ]

    mass_label = tasks.mass
    mass = str(submit_config.mass[mass_label])
    modules.append(templates.action(name=f"stag_mass_{mass_label}",
                                    mass=mass,
                                    gauge_fat='gauge_fat',
                                    gauge_long='gauge_long'))

    nvecs = str(3 * submit_config.time)
    w_indices = tasks.w_indices if tasks.w_indices else list(range(submit_config.noise))
    v_indices = tasks.v_indices if tasks.v_indices else list(range(submit_config.noise))
    pairings = list(sorted(filter(lambda x: x[0] != x[1], set((tuple(x) for x in itertools.product(w_indices,v_indices))))))
    module_set = set()
    vec_path = outfile_config_list.a2avec.filestem
    meson_path = outfile_config_list.meson.filestem

    for w_index, v_index in pairings:
        v_name = f"v{v_index}"
        w_name = f"w{w_index}"

        if w_name not in module_set:
            module_set.add(w_name)
            modules.append(templates.time_diluted_noise(w_name, 1))

        if not tasks.low_memory_mode:
            v_name_unique = v_name
        else:
            v_name_unique = v_name + f"_{w_name}"

        if v_name_unique not in module_set:
            module_set.add(v_name_unique)
            infile = vec_path.format(mass=submit_config.mass_out_label[mass_label],
            seed_index=str(v_index),
            **submit_conf_dict)

            modules.append(templates.load_vectors(name=v_name_unique,
                                                  filestem=infile,
                                                  size=nvecs,
                                                  multifile='true' if tasks.multifile else 'false'))

        outfile = meson_path.format(mass=submit_config.mass_out_label[mass_label],
                                    w_index=w_index,
                                    v_index=v_index,
                                    **submit_conf_dict)

        modules.append(templates.meson_field(name=f'mf_{w_index}_{v_index}',
                                             action=f"stag_mass_{mass_label}",
                                             block=nvecs,
                                             gammas=tasks.gammas,
                                             gauge='',
                                             low_modes='',
                                             left=w_name+"_vec",
                                             right=v_name_unique,
                                             output=outfile,
                                             apply_g5='false'))

    schedule = [m["id"]['name'] for m in modules]

    return modules, schedule


def bad_files(submit_config: SubmitHadronsConfig,
              task_config: TaskBase, outfile_config_list: OutfileList) -> t.List[str]:
    logging.warning(
        "Check completion succeeds automatically. No implementation of bad_files function in `hadrons_a2a_vectors.py`.")
    return []


def get_task_factory():
    return A2ASIBTask
