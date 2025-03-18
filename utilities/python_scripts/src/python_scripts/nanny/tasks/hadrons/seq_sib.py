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
import logging
import typing as t

from pydantic.dataclasses import dataclass

from python_scripts import Gamma
from python_scripts.nanny import TaskBase
from python_scripts.nanny.config import OutfileList
from python_scripts.nanny.tasks.hadrons import SubmitHadronsConfig, templates


@dataclass
class SeqSIBTask(TaskBase):
    mass: str
    gammas: t.List[Gamma]
    free: bool = False

    # def __init__(self,mass: str, gammas: t.List[str]):
    #     self.mass = mass
    #     self.gammas = [Gamma[g] for g in gammas]

def input_params(tasks: SeqSIBTask, submit_config: SubmitHadronsConfig, outfile_config_list: OutfileList) -> t.Tuple[t.List[t.Dict], t.Optional[t.List[str]]]:

    submit_conf_dict = submit_config.string_dict()

    run_tsources = list(map(str, submit_config.tsource_range))

    gauge_filepath = outfile_config_list.gauge_links.filestem.format(**submit_conf_dict)
    gauge_fat_filepath = outfile_config_list.fat_links.filestem.format(**submit_conf_dict)
    gauge_long_filepath = outfile_config_list.long_links.filestem.format(**submit_conf_dict)

    if tasks.free:
        modules = [
            templates.unit_gauge('gauge'),
            templates.unit_gauge('gauge_fat'),
            templates.unit_gauge('gauge_long'),
            templates.cast_gauge('gauge_fatf', 'gauge_fat'),
            templates.cast_gauge('gauge_longf', 'gauge_long')
        ]
    else:
        gauge_filepath = outfile_config_list.gauge_links.filestem.format(**submit_conf_dict)
        gauge_fat_filepath = outfile_config_list.fat_links.filestem.format(**submit_conf_dict)
        gauge_long_filepath = outfile_config_list.long_links.filestem.format(**submit_conf_dict)

        modules = [
            templates.load_gauge('gauge', gauge_filepath),
            templates.load_gauge('gauge_fat', gauge_fat_filepath),
            templates.load_gauge('gauge_long', gauge_long_filepath),
            templates.cast_gauge('gauge_fatf', 'gauge_fat'),
            templates.cast_gauge('gauge_longf', 'gauge_long')
        ]

    mass_label = tasks.mass
    name = f"stag_mass_{mass_label}"
    mass = str(submit_config.mass[mass_label])
    modules.append(templates.action(name=name,
                                    mass=mass,
                                    gauge_fat='gauge_fat',
                                    gauge_long='gauge_long'))

    name = f"istag_mass_{mass_label}"
    mass = str(submit_config.mass[mass_label])
    modules.append(templates.action_float(name=name,
                                          mass=mass,
                                          gauge_fat='gauge_fatf',
                                          gauge_long='gauge_longf'))

    multifile = 'false'
    epack_path = outfile_config_list.eig.filestem.format(**submit_conf_dict)

    modules.append(templates.epack_load(name='epack',
                                        filestem=epack_path,
                                        size=submit_conf_dict['sourceeigs'],
                                        multifile=multifile))

    mass = str(submit_config.mass[mass_label])
    modules.append(templates.epack_modify(name=f"evecs_mass_{mass_label}",
                                          eigen_pack='epack',
                                          mass=mass))

    modules.append(templates.mixed_precision_cg(
        name=f"stag_ama_mass_{mass_label}",
        outer_action=f"stag_mass_{mass_label}",
        inner_action=f"istag_mass_{mass_label}",
        residual=submit_conf_dict['cg_residual']
    ))

    modules.append(templates.lma_solver(
        name=f"stag_ranLL_mass_{mass_label}",
        action=f"stag_mass_{mass_label}",
        low_modes=f"evecs_mass_{mass_label}"
    ))
    modules.append(templates.sink(name='sink', mom='0 0 0'))

    for tsource in run_tsources:
        modules.append(templates.noise_rw(
            name=f"noise_t{tsource}",
            nsrc=submit_conf_dict['noise'],
            t0=tsource,
            tstep=submit_conf_dict['time']
        ))

        for i, gamma in enumerate(tasks.gammas):
            glabel = gamma.name

            if i == 0:
                assert gamma == Gamma.G5_G5

            high_path = outfile_config_list.seq_modes.filestem
            solver_labels = ['ranLL', 'ama']
            for slabel in solver_labels:

                quark = f"quark_{slabel}_mass_{mass_label}_t{tsource}_{glabel}"

                source = f"noise_t{tsource}"
                solver = f"stag_{slabel}_mass_{mass_label}"
                if slabel == "ama":
                    guess = f"quark_ranLL_mass_{mass_label}_t{tsource}_{glabel}"
                else:
                    guess = ""

                modules.append(templates.quark_prop(
                    name=quark,
                    source=source,
                    solver=solver,
                    guess=guess,
                    gammas=gamma.gamma_string,
                    apply_g5='true',
                    gauge="" if gamma.local else "gauge"
                ))

                # sib = Gamma.G1_G1
                # quark_g1 = quark+sib.name
                # modules.append(templates.seq_gamma(
                #     name=quark_g1,
                #     q=quark+glabel,
                #     ta='0',
                #     tb=submit_conf_dict['time'],
                #     gammas=sib.gamma_string,
                #     apply_g5='false',
                #     gauge="",
                #     mom='0 0 0'
                # ))

                quark_g1_M = quark + "_M"
                guess_M = guess + "_M" if guess else guess
                modules.append(templates.quark_prop(
                    name=quark + "_M",
                    source=quark+glabel,
                    solver=solver,
                    guess=guess_M,
                    gammas="",
                    apply_g5='false',
                    gauge=""
                ))

                pion = f"quark_{slabel}_mass_{mass_label}_t{tsource}_G5_G5G5_G5"

                output = high_path.format(
                    mass=submit_config.mass_out_label[mass_label],
                    dset=slabel,
                    gamma=glabel,
                    tsource=tsource,
                    **submit_conf_dict
                )

                modules.append(templates.prop_contract(
                    name=f"corr_{slabel}_{glabel}_mass_{mass_label}_t{tsource}",
                    source=quark_g1_M,
                    sink=pion,
                    sink_func='sink',
                    source_shift=f"noise_t{tsource}_shift",
                    source_gammas="",
                    sink_gammas=gamma.gamma_string,
                    apply_g5='true',
                    gauge="" if gamma.local else "gauge",
                    output=output
                ))

    schedule = [m["id"]['name'] for m in modules]

    return modules, schedule


def bad_files(task_config: TaskBase, submit_config: SubmitHadronsConfig,
              outfile_config_list: OutfileList) -> t.List[str]:
    logging.warning(
        "Check completion succeeds automatically. No implementation of bad_files function in `hadrons_a2a_vectors.py`.")
    return []


def get_task_factory():
    return SeqSIBTask.from_dict
