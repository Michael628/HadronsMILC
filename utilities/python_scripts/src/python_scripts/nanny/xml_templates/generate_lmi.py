import itertools
import re

from python_scripts.utils import deep_copy_dict
from python_scripts.nanny.xml_templates.components import (
    gauge, eig, meson, highmodes, wrapper
)
import typing as t


GAMMAS = ['pion_local', 'vec_local', 'vec_onelink']


def append_gauges_and_actions(modules: t.List, params: dict) -> None:

    modules += gauge.gauge_dp_params(
        ens=params['ens'],
        series=params['series'],
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

    for mass_label, mass in params['mass'].items():

        modules += gauge.action_params(
            name=f"stag_mass_{mass_label}",
            mass=mass,
            fat_label='gauge_fat',
            long_label='gauge_long',
            single_precision=False
        )

        modules += gauge.action_params(
            name=f"istag_mass_{mass_label}",
            mass=mass,
            fat_label='gauge_fatf',
            long_label='gauge_longf',
            single_precision=True
        )


def append_epacks(modules: t.List,
                  load: bool,
                  params: dict,
                  save_evals: bool = False) -> None:

    modules += gauge.op_params('stag_op', 'stag_e')

    if load:
        modules += eig.epack_load_params(
            name='epack',
            filestem=(f"eigen/eig{params['ens']}nv"
                      f"{params['sourceeigs']}{params['series']}"),
            eigs=params['eigs'],
            multifile=params['multifile']
        )
    else:
        modules += eig.irl_params(
            name='epack',
            operator='stag_op_schur',
            lanczos_params={
                'alpha': params['alpha'],
                'beta': params['beta'],
                'npoly': params['npoly'],
                'nstop': params['nstop'],
                'nk': params['nk'],
                'nm': params['nm']
            }
        )

    for mass_label, mass in params['mass'].items():

        modules += eig.epack_modify_params(
            name=f"evecs_mass_{mass_label}",
            epack="epack",
            mass=mass
        )

    if save_evals:
        modules += eig.eval_save_params(
            output=(f"eigen/eval/eval{params['ens']}nv"
                    f"{params['sourceeigs']}{params['series']}")
        )


def append_mesons(modules: t.List,
                  operators: t.List[str], params: dict) -> None:

    if 'l' not in params['mass']:
        raise ValueError(("Masses must include light "
                          "quark key `l` for meson generation"))

    mass_label = "mass_l"
    mass_output = f"m{str(params['mass']['l'])[2:]}"

    if 'local' in operators:
        modules += meson.meson_local_params(
            name=f"mf_local_{mass_label}",
            action=f"stag_{mass_label}",
            block=params['blocksize'],
            lowmodes=f"evecs_{mass_label}",
            output=(f"e{params['eigs']}n{params['noise']}dt{params['dt']}/"
                    f"mesons/{mass_output}/mf_{params['series']}"),
            include_pion=True
        )

    if 'onelink' in operators:
        modules += meson.meson_onelink_params(
            name=f"mf_onelink_{mass_label}",
            action=f"stag_{mass_label}",
            block=params['blocksize'],
            lowmodes=f"evecs_{mass_label}",
            gauge='gauge',
            output=(f"e{params['eigs']}n{params['noise']}dt{params['dt']}/"
                    f"mesons/{mass_output}/mf_{params['series']}")
        )


def append_solvers(modules: t.List, eigs: bool, params: dict) -> None:

    for mass_label in params['mass'].keys():

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


def append_sources_and_quarks(modules: t.List, eigs: bool, params: dict) \
        -> None:

    for tslice in params['time_range']:
        modules += highmodes.source_params(
            name=f"noise_t{tslice}",
            nsrc=params['noise'],
            tstep=params['time'],
            t0=str(tslice)
        )

    def m1_eq_m2(x):
        return x[-2] == x[-1]

    for (tslice, glabel, slabel, mlabel, _) in filter(m1_eq_m2,
                                                      params['quark_iter']):

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


def append_contractors(modules: t.List, params: dict) -> None:

    modules += highmodes.sink_params()

    def m1_ge_m2(x):
        return x[-2] >= x[-1]

    for (tslice, glabel, slabel, m1label, m2label) \
            in filter(m1_ge_m2, params['quark_iter']):

        quark1 = f"quark_{slabel}_{glabel}_mass_{m1label}_t{tslice}"
        quark2 = f"quark_{slabel}_pion_local_mass_{m2label}_t{tslice}"

        if m1label == m2label:
            mass_label = f"mass_{m1label}"
            mass_output = f"m{str(params['mass'][m1label])[2:]}"
        else:
            mass_label = f"mass_{m1label}_mass_{m2label}"
            mass_output = (f"m{str(params['mass'][m1label])[2:]}"
                           f"_m{str(params['mass'][m2label])[2:]}")

        module_builder = getattr(highmodes, f"contractor_{glabel}_params")
        kwargs = {
            'name': f"corr_{slabel}_{glabel}_{mass_label}_t{tslice}",
            'source': quark1,
            'sink': quark2,
            'source_shift': f"noise_t{tslice}_shift",
            'output': (f"e{params['eigs']}n{params['noise']}dt{params['dt']}/"
                       f"correlators/{mass_output}/{glabel}/{slabel}/"
                       f"corr_{glabel}_{slabel}_{mass_output}_t{tslice}"
                       f"_{params['series']}")
        }
        if 'local' not in glabel:
            kwargs['gauge'] = "gauge"

        modules += module_builder(**kwargs)


def build_schedule(module_info, params):

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
        for i, mass in enumerate(params['mass'].keys()):
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


def build_params(tasks: t.Dict, **kwargs):

    modules = []

    params = deep_copy_dict(kwargs)

    xml_params = wrapper(
        runid=(f"LMI-RW-series-{params['series']}"
               f"-{params['eigs']}-eigs-{params['noise']}-noise"),
        sched=params.get('schedule', ''),
        cfg=params['cfg']
    )

    epack = tasks.get('epack', {})
    has_eigs = len(epack) != 0

    t_start = int(params["tstart"])
    t_stop = int(params["tstop"]) + 1
    t_step = int(params["dt"])
    params['time_range'] = list(range(t_start, t_stop, t_step))

    append_gauges_and_actions(modules, params)

    if has_eigs:
        append_epacks(modules, params=params, **epack)

    meson_operators = tasks.get('meson', [])
    append_mesons(modules, meson_operators, params)

    high_params = tasks.get('high_modes', {})
    if high_params:
        assert all([g in GAMMAS for g in high_params.keys()])
        assert 'pion_local' in high_params.keys()
        solve_combos = []
        solver_labels = ["ranLL", "ama"] if epack else ['ama']

        for gamma, solve_params in high_params.items():
            if 'mass' not in solve_params:
                raise ValueError(("must include `mass` param"
                                  " for high mode solves"))

            mass_labels: t.Union[str, t.List[str]] = solve_params['mass']
            if isinstance(mass_labels, str):
                if mass_labels.strip().lower() == 'all':
                    mass_labels = list(params['mass'].keys())
                else:
                    mass_labels = [mass_labels]

            assert all([m in params['mass'].keys() for m in mass_labels])

            solve_combos.append(itertools.product(
                params['time_range'],
                [gamma],
                solver_labels,
                mass_labels,
                mass_labels
            ))

        params['quark_iter'] = list(itertools.chain(*solve_combos))

        append_solvers(modules, has_eigs, params)
        append_sources_and_quarks(modules, has_eigs, params)
        append_contractors(modules, params)

    xml_params['grid']["modules"] = {"module": modules}

    module_info = [m["id"] for m in modules]
    schedule = build_schedule(module_info, params)

    return xml_params, schedule

    # return xml_params
