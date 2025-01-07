import itertools
import re

from python_scripts.nanny import config
from python_scripts.nanny.xml_templates import gauge, meson, highmodes, eig
import typing as t


GAMMAS = ['pion_local', 'vec_local', 'vec_onelink']


def append_gauges_and_actions(modules: t.List, run_config: config.RunConfig) -> None:

    modules += gauge.gauge_dp_params(
        ens=run_config.ens,
        series=run_config.series,
        include_unsmeared=True
    )

    modules += gauge.gauge_sp_params(
        field_names=["gauge_fat", "gauge_long"]
    )

    modules += gauge.action_params(
        name='stag_e',
        mass='0.0',
        fat_label='gauge_fat',
        long_label='gauge_long'
    )

    for mass_label, mass in run_config.mass.items():

        modules += gauge.action_params(
            name=f"stag_mass_{mass_label}",
            mass=str(mass),
            fat_label='gauge_fat',
            long_label='gauge_long',
            single_precision=False
        )

        modules += gauge.action_params(
            name=f"istag_mass_{mass_label}",
            mass=str(mass),
            fat_label='gauge_fatf',
            long_label='gauge_longf',
            single_precision=True
        )


def append_epacks(modules: t.List,
                  load: bool,
                  run_config: config.RunConfig,
                  save_evals: bool = False,
                  save_eigs: bool = False) -> None:

    modules += gauge.op_params('stag_op', 'stag_e')

    if load:
        modules += eig.epack_load_params(
            name='epack',
            filestem=(f"eigen/eig{run_config.ens}nv"
                      f"{run_config.sourceeigs}{run_config.series}"),
            eigs=str(run_config.eigs),
            multifile=str(run_config.multifile)
        )
    else:
        eig_out: str = ""
        if save_eigs:
            eig_out = (f"eig{run_config.ens}nv{run_config.sourceeigs}"
                       f"{run_config.series}")
        modules += eig.irl_params(
            name='epack',
            operator='stag_op_schur',
            lanczos_params={
                'alpha': run_config.alpha,
                'beta': run_config.beta,
                'npoly': run_config.npoly,
                'nstop': run_config.nstop,
                'nk': run_config.nk,
                'nm': run_config.nm
            }
        )

    for mass_label, mass in run_config.mass.items():

        modules += eig.epack_modify_params(
            name=f"evecs_mass_{mass_label}",
            epack="epack",
            mass=str(mass)
        )

    if save_evals:
        modules += eig.eval_save_params(
            output=(f"eigen/eval/eval{run_config.ens}nv"
                    f"{run_config.sourceeigs}{run_config.series}")
        )


def append_mesons(modules: t.List,
                  operators: t.List[str], run_config: config.RunConfig) -> None:

    if 'l' not in run_config.mass:
        raise ValueError(("Masses must include light "
                          "quark key `l` for meson generation"))

    mass_label = "mass_l"
    mass_output = f"m{str(run_config.mass['l'])[2:]}"

    if 'local' in operators:
        modules += meson.meson_local_params(
            name=f"mf_local_{mass_label}",
            action=f"stag_{mass_label}",
            block=str(run_config.blocksize),
            lowmodes=f"evecs_{mass_label}",
            output=(f"e{run_config.eigs}n{run_config.noise}dt{run_config.dt}/"
                    f"mesons/{mass_output}/mf_{run_config.series}"),
            include_pion=True
        )

    if 'onelink' in operators:
        modules += meson.meson_onelink_params(
            name=f"mf_onelink_{mass_label}",
            action=f"stag_{mass_label}",
            block=str(run_config.blocksize),
            lowmodes=f"evecs_{mass_label}",
            gauge='gauge',
            output=(f"e{run_config.eigs}n{run_config.noise}dt{run_config.dt}/"
                    f"mesons/{mass_output}/mf_{run_config.series}")
        )


def append_solvers(modules: t.List, eigs: bool, run_config: config.RunConfig) -> None:

    for mass_label in run_config.mass.keys():

        modules += highmodes.mixed_precision_solver_params(
            name=f"stag_ama_mass_{mass_label}",
            outer_action=f"stag_mass_{mass_label}",
            inner_action=f"istag_mass_{mass_label}",
            residual='1e-8'
        )

        if eigs:
            modules += highmodes.ranLL_solver_params(
                name=f"stag_ranLL_mass_{mass_label}",
                action=f"stag_mass_{mass_label}",
                lowmodes=f"evecs_mass_{mass_label}"
            )


def append_sources_and_quarks(modules: t.List, eigs: bool, quark_iter: t.List, run_config: config.RunConfig) \
        -> None:

    for tslice in run_config.time_range:
        modules += highmodes.source_params(
            name=f"noise_t{tslice}",
            nsrc=str(run_config.noise),
            tstep=str(run_config.time),
            t0=str(tslice)
        )

    def m1_eq_m2(x):
        return x[-2] == x[-1]

    for (tslice, glabel, slabel, mlabel, _) in filter(m1_eq_m2,
                                                      quark_iter):

        quark = f"quark_{slabel}_{glabel}_mass_{mlabel}_t{tslice}"
        source = f"noise_t{tslice}"
        solver = f"stag_{slabel}_mass_{mlabel}"
        if slabel == "ama" and eigs:
            guess = f"quark_ranLL_{glabel}_mass_{mlabel}_t{tslice}"
        else:
            guess = ""

        module_builder = getattr(highmodes, f"propagator_{glabel}_params")
        kwargs = {
            'name': quark,
            'source': source,
            'solver': solver,
            'guess': guess
        }
        if 'local' not in glabel:
            kwargs['gauge'] = "gauge"

        modules += module_builder(**kwargs)


def append_contractors(modules: t.List, quark_iter: t.List, run_config: config.RunConfig) -> None:

    modules += highmodes.sink_params()

    def m1_ge_m2(x):
        return x[-2] >= x[-1]

    for (tslice, glabel, slabel, m1label, m2label) \
            in filter(m1_ge_m2, quark_iter):

        quark1 = f"quark_{slabel}_{glabel}_mass_{m1label}_t{tslice}"
        quark2 = f"quark_{slabel}_pion_local_mass_{m2label}_t{tslice}"

        if m1label == m2label:
            mass_label = f"mass_{m1label}"
            mass_output = f"m{str(run_config.mass[m1label])[2:]}"
        else:
            mass_label = f"mass_{m1label}_mass_{m2label}"
            mass_output = (f"m{str(run_config.mass[m1label])[2:]}"
                           f"_m{str(run_config.mass[m2label])[2:]}")

        module_builder = getattr(highmodes, f"contractor_{glabel}_params")
        kwargs = {
            'name': f"corr_{slabel}_{glabel}_{mass_label}_t{tslice}",
            'source': quark1,
            'sink': quark2,
            'source_shift': f"noise_t{tslice}_shift",
            'output': (f"e{run_config.eigs}n{run_config.noise}dt{run_config.dt}/"
                       f"correlators/{mass_output}/{glabel}/{slabel}/"
                       f"corr_{glabel}_{slabel}_{mass_output}_t{tslice}"
                       f"_{run_config.series}")
        }
        if 'local' not in glabel:
            kwargs['gauge'] = "gauge"

        modules += module_builder(**kwargs)


def build_schedule(module_info, run_config):

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
        lambda x: "mass" not in x['name'] and "_t" not in x['name']
    )

    sorted_modules = dp_gauges + indep_mass_tslice
    sorted_modules += meson_field_inputs + meson_fields
    sorted_modules += sp_gauges

    def gamma_order(x):
        for i, gamma in enumerate(GAMMAS):
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


def build_xml_params(tasks: t.Dict, run_config: config.RunConfig):

    modules = []

    epack = tasks.get('epack', {})
    has_eigs = len(epack) != 0

    append_gauges_and_actions(modules, run_config)

    if has_eigs:
        append_epacks(modules, run_config=run_config, **epack)

    meson_operators = tasks.get('meson', [])
    append_mesons(modules, meson_operators, run_config)

    high_params = tasks.get('high_modes', {})
    if high_params:
        high_config = config.get_config_factory('high_modes')(high_params)

        solve_combos = []
        solver_labels = ["ranLL", "ama"] if epack else ['ama']

        for op in high_config.gammas:
            mass_labels: t.List[str] = op.mass

            assert all([m in run_config.mass.keys() for m in mass_labels])

            solve_combos.append(itertools.product(
                run_config.time_range,
                [op.gamma.name.lower()],
                solver_labels,
                mass_labels,
                mass_labels
            ))

        quark_iter = list(itertools.chain(*solve_combos))

        append_solvers(modules, has_eigs, run_config)
        append_sources_and_quarks(modules, has_eigs, quark_iter, run_config)
        append_contractors(modules, quark_iter, run_config)

    module_info = [m["id"] for m in modules]
    schedule = build_schedule(module_info, run_config)

    return modules, schedule

    # return xml_params
