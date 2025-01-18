import itertools
import os.path
import re

import pandas as pd

from python_scripts.nanny.runio import hadrons
from python_scripts import Gamma, utils
from python_scripts.nanny.config import OutfileList, SubmitHadronsConfig, TaskBase
import typing as t

def build_params(submit_config: SubmitHadronsConfig, tasks: TaskBase,
                 outfile_config_list: OutfileList) -> t.Tuple[t.List[t.Dict], t.Optional[t.List[str]]]:

     modules = []

     submit_conf_dict = submit_config.string_dict()

     run_tsources = list(map(str, submit_config.tsource_range))

     gauge_filepath = outfile_config_list.gauge_links.filestem.format(**submit_conf_dict)
     gauge_fat_filepath = outfile_config_list.fat_links.filestem.format(**submit_conf_dict)
     gauge_long_filepath = outfile_config_list.long_links.filestem.format(**submit_conf_dict)

     modules = [
          hadrons.load_gauge('gauge', gauge_filepath),
          hadrons.load_gauge('gauge_fat', gauge_fat_filepath),
          hadrons.load_gauge('gauge_long', gauge_long_filepath),
          hadrons.cast_gauge('gauge_fatf', 'gauge_fat'),
          hadrons.cast_gauge('gauge_longf', 'gauge_long')
     ]

     epack_path = outfile_config_list.eig.filestem.format(**submit_conf_dict)

     # Load or generate eigenvectors
     modules.append(hadrons.epack_load(name='epack',
                                       filestem=epack_path,
                                       size=submit_conf_dict['sourceeigs'],
                                       multifile=submit_conf_dict['multifile']))

     mass_label = 'l'
     name = f"stag_mass_{mass_label}"
     mass = str(submit_config.mass[mass_label])
     modules.append(hadrons.action(name=name,
                                   mass=mass,
                                   gauge_fat='gauge_fat',
                                   gauge_long='gauge_long'))

     modules.append(hadrons.epack_modify(name=f"evecs_mass_{mass_label}",
                                         eigen_pack='epack',
                                         mass=mass))

     name = f"istag_mass_{mass_label}"
     modules.append(hadrons.action_float(name=name,
                                         mass=mass,
                                         gauge_fat='gauge_fatf',
                                         gauge_long='gauge_longf'))

     modules.append(hadrons.lma_solver(
          name=f"stag_ranLL_mass_{mass_label}",
          action=f"stag_mass_{mass_label}",
          low_modes=f"evecs_mass_{mass_label}"
     ))

     modules.append(hadrons.mixed_precision_cg(
          name=f"stag_ama_mass_{mass_label}",
          outer_action=f"stag_mass_{mass_label}",
          inner_action=f"istag_mass_{mass_label}",
          residual='1e-8'
     ))

     gamma_label="G1_G1"
     gamma="(G1 G1)"
     modules.append(hadrons.a2a_vector(
          name="a2avec",
          noise='noise',
          action=f'stag_{mass_label}',
          low_modes=f"quark_ranLL_{gamma_label}_{mass_string}" + gamma_label,
          solver,
          high_output
     ))

     module_info = [m["id"] for m in modules]
     schedule = build_schedule(module_info)

     module = copy.deepcopy(module_templates["a2aVector"])
     module["options"]["action"] = f"stag_{mass_string}"
     module["options"]["lowModes"] =
     module["options"]["solver"] = f"stag_ama_{mass_string}"
     module["options"]["highOutput"] = f"data/vectors/{env['SEED']}"
     module["options"]["highMultiFile"] = "true"
     module["options"]["norm2"] = env['NOISE']
     modules.append(module)


     params["grid"]["modules"] = {"module":modules}
     
     moduleList = [m["id"]["name"] for m in modules]

     f = open(schedule_file, "w")
     f.write(str(len(moduleList)) + "\n" + "\n".join(moduleList))
     f.close()

     return params
