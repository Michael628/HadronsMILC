import itertools
import re

from python_scripts.nanny import xml_templates
from python_scripts import ConfigBase, Gamma
from python_scripts.nanny.config import OutfileConfigList, RunConfig, GenerateLMITaskConfig
import typing as t


def build_schedule(module_info, run_config):
    gammas = ['pion_local', 'vec_local', 'vec_onelink']

    def pop_conditional(mi, cond):
        indices = [
            i
            for i, item in enumerate(mi)
            if cond(item)
        ]
        # Pop in reverse order but return in original order
        return [mi.pop(i) for i in indices[::-1]][::-1]

    def get_mf_inputs(x):
        is_action = x['type'].endswith('ImprovedStaggeredMILC')
        is_evec = 'ModifyEigenPack' in x['type']
        has_sea_mass = "mass_l" in x['name']
        return (is_action or is_evec) and has_sea_mass

    dp_gauges = pop_conditional(
        module_info,
        lambda x: 'LoadIldg' in x['type']
    )
    sp_gauges = pop_conditional(
        module_info,
        lambda x: 'PrecisionCast' in x['type']
    )
    meson_fields = pop_conditional(
        module_info,
        lambda x: 'A2AMesonField' in x['type']
    )
    meson_field_inputs = pop_conditional(module_info, get_mf_inputs)

    indep_mass_tslice = pop_conditional(
        module_info,
        lambda x: ("mass" not in x['name'] or "mass_zero" in x['name']) and "_t" not in x['name']
    )

    sorted_modules = dp_gauges + indep_mass_tslice
    sorted_modules += meson_field_inputs + meson_fields
    sorted_modules += sp_gauges

    def gamma_order(x):
        for i, gamma in enumerate(gammas):
            if gamma in x['name']:
                return i
        return -1

    def mass_order(x):
        for i, mass in enumerate(run_config.mass.keys()):
            if f"mass_{mass}" in x['name']:
                return i
        return -1

    def mixed_mass_last(x):
        return len(re.findall(r'_mass', x['name']))

    def tslice_order(x):
        time = re.findall(r'_t(\d+)', x['name'])
        if len(time):
            return int(time[0])
        else:
            return -1

    module_info = sorted(module_info, key=gamma_order)
    module_info = sorted(module_info, key=mass_order)
    module_info = sorted(module_info, key=mixed_mass_last)
    module_info = sorted(module_info, key=tslice_order)

    sorted_modules += module_info

    return [m['name'] for m in sorted_modules]


def build_xml_params(tasks: GenerateLMITaskConfig, run_config: RunConfig, outfile_config: OutfileConfigList):
    run_conf_dict = run_config.string_dict

    gauge_filepath = outfile_config.gauge_links.filestem.format(**run_conf_dict)
    gauge_fat_filepath = outfile_config.fat_links.filestem.format(**run_conf_dict)
    gauge_long_filepath = outfile_config.long_links.filestem.format(**run_conf_dict)

    modules = [
        xml_templates.load_gauge('gauge', gauge_filepath),
        xml_templates.load_gauge('gauge_fat', gauge_fat_filepath),
        xml_templates.load_gauge('gauge_long', gauge_long_filepath),
        xml_templates.cast_gauge('gauge_fatf', 'gauge_fat'),
        xml_templates.cast_gauge('gauge_longf', 'gauge_long')
    ]

    for mass_label in tasks.mass:
        name = f"stag_mass_{mass_label}"
        mass = str(run_config.mass[mass_label])
        modules.append(xml_templates.action(name=name,
                                            mass=mass,
                                            gauge_fat='gauge_fat',
                                            gauge_long='gauge_long'))

    if tasks.high_modes:
        for mass_label in tasks.high_modes.mass:
            name = f"istag_mass_{mass_label}"
            mass = str(run_config.mass[mass_label])
            modules.append(xml_templates.action_float(name=name,
                                                      mass=mass,
                                                      gauge_fat='gauge_fatf',
                                                      gauge_long='gauge_longf'))

    if tasks.epack:
        epack_path = ''
        if tasks.epack.load or tasks.epack.save_eigs:
            epack_path = outfile_config.eig.filestem.format(**run_conf_dict)

        if tasks.epack.load:
            modules.append(xml_templates.epack_load(name='epack',
                                                    filestem=epack_path,
                                                    size=run_conf_dict['sourceeigs'],
                                                    multifile=run_conf_dict['multifile']))
        else:
            modules.append(xml_templates.op('stag_op', 'stag_mass_zero'))
            modules.append(xml_templates.irl(name='epack',
                                             op='stag_op_schur',
                                             alpha=run_conf_dict['alpha'],
                                             beta=run_conf_dict['beta'],
                                             npoly=run_conf_dict['npoly'],
                                             nstop=run_conf_dict['nstop'],
                                             nk=run_conf_dict['nk'],
                                             nm=run_conf_dict['nm'],
                                             output=epack_path))

        for mass_label in tasks.mass:
            if mass_label == "zero":
                continue
            mass = str(run_config.mass[mass_label])
            modules.append(xml_templates.epack_modify(name=f"evecs_mass_{mass_label}",
                                                      eigen_pack='epack',
                                                      mass=mass))

        if tasks.epack.save_evals:
            eval_path = outfile_config.eval.filestem.format(**run_conf_dict)
            modules.append(xml_templates.eval_save(name='eval_save',
                                                   eigen_pack='epack',
                                                   output=eval_path))

    if tasks.meson:
        meson_path = outfile_config.meson.filestem
        for op in tasks.meson.operations:
            op_type = op.gamma.name.lower()
            gauge = "" if op.gamma == Gamma.LOCAL else "gauge"
            for mass_label in op.mass:
                output = meson_path.format(
                    mass=run_config.mass_out_label[mass_label],
                    **run_conf_dict
                )
                modules.append(xml_templates.meson_field(
                    name=f"mf_{op_type}_mass_{mass_label}",
                    action=f"stag_mass_{mass_label}",
                    block=run_conf_dict['blocksize'],
                    gammas=op.gamma.gamma_string,
                    apply_g5='false',
                    gauge=gauge,
                    low_modes=f"evecs_mass_{mass_label}",
                    left="",
                    right="",
                    output=output
                ))

    if tasks.high_modes:

        modules.append(xml_templates.sink(name='sink', mom='0 0 0'))

        for mass_label in tasks.high_modes.mass:
            modules.append(xml_templates.mixed_precision_cg(
                name=f"stag_ama_mass_{mass_label}",
                outer_action=f"stag_mass_{mass_label}",
                inner_action=f"istag_mass_{mass_label}",
                residual='1e-8'
            ))

            if tasks.epack:
                modules.append(xml_templates.lma_solver(
                    name=f"stag_ranLL_mass_{mass_label}",
                    action=f"stag_mass_{mass_label}",
                    low_modes=f"evecs_mass_{mass_label}"
                ))

        for tslice in map(str, run_config.time_range):
            modules.append(xml_templates.noise_rw(
                name=f"noise_t{tslice}",
                nsrc=run_conf_dict['noise'],
                t0=tslice,
                tstep=run_conf_dict['time']
            ))

        def m1_eq_m2(x):
            return x[-2] == x[-1]

        solver_labels = ["ranLL", "ama"] if tasks.epack else ['ama']

        high_path = outfile_config.high_modes.filestem
        for op in tasks.high_modes.operations:
            glabel = op.gamma.name.lower()
            quark_iter = list(itertools.product(
                list(map(str, run_config.time_range)),
                solver_labels,
                op.mass,
                op.mass
            ))

            for (tslice, slabel, mlabel, _) in filter(m1_eq_m2, quark_iter):

                quark = f"quark_{slabel}_{glabel}_mass_{mlabel}_t{tslice}"
                source = f"noise_t{tslice}"
                solver = f"stag_{slabel}_mass_{mlabel}"
                if slabel == "ama" and tasks.epack:
                    guess = f"quark_ranLL_{glabel}_mass_{mlabel}_t{tslice}"
                else:
                    guess = ""

                modules.append(xml_templates.quark_prop(
                    name=quark,
                    source=source,
                    solver=solver,
                    guess=guess,
                    gammas=op.gamma.gamma_string,
                    apply_g5='true',
                    gauge="" if op.gamma.local else "gauge"
                ))

            def m1_ge_m2(x):
                return x[-2] >= x[-1]

            for (tslice, slabel, m1label, m2label) \
                    in filter(m1_ge_m2, quark_iter):

                quark1 = f"quark_{slabel}_{glabel}_mass_{m1label}_t{tslice}"
                quark2 = f"quark_{slabel}_pion_local_mass_{m2label}_t{tslice}"

                if m1label == m2label:
                    mass_label = f"mass_{m1label}"
                    mass_output = f"{run_config.mass_out_label[m1label]}"
                else:
                    mass_label = f"mass_{m1label}_mass_{m2label}"
                    mass_output = (f"{run_config.mass_out_label[m1label]}"
                                   f"_m{run_config.mass_out_label[m2label]}")

                output = high_path.format(
                    mass=mass_output,
                    dset=slabel,
                    gamma=glabel,
                    tsource=tslice,
                    **run_conf_dict
                )

                modules.append(xml_templates.prop_contract(
                    name=f"corr_{slabel}_{glabel}_{mass_label}_t{tslice}",
                    source=quark1,
                    sink=quark2,
                    sink_func='sink',
                    source_shift=f"noise_t{tslice}_shift",
                    source_gammas=op.gamma.gamma_string,
                    sink_gammas=op.gamma.gamma_string,
                    apply_g5='true',
                    gauge="" if op.gamma.local else "gauge",
                    output=output
                ))

    module_info = [m["id"] for m in modules]
    schedule = build_schedule(module_info, run_config)

    return modules, schedule


def generate_outfile_formatting(task_config: ConfigBase, outfile_config: OutfileConfigList, run_config: RunConfig):
    assert isinstance(task_config, GenerateLMITaskConfig)

    if task_config.epack:
        if task_config.epack.save_eigs:
            if task_config.epack.multifile:
                yield {'eig_index': list(range(int(run_config.eigs)))}, outfile_config.eigdir
            else:
                yield {}, outfile_config.eig
        if task_config.epack.save_eigs:
            yield {}, outfile_config.eval

    if task_config.meson:
        res: t.Dict = {}
        for op in task_config.meson.operations:
            res['gamma'] = op.gamma.gamma_list
            res['mass'] = [run_config.mass_out_label[m] for m in op.mass]
            yield res, outfile_config.meson

    if task_config.high_modes:
        res = {'tsource': list(range(run_config.tstart, run_config.tstop, run_config.dt))}
        if task_config.epack:
            res['dset'] = ['ama', 'ranLL']
        else:
            res['dset'] = ['ama']

        for op in task_config.high_modes.operations:
            res['gamma'] = op.gamma.name.lower()
            res['mass'] = [run_config.mass_out_label[m] for m in op.mass]
            yield res, outfile_config.high_modes
