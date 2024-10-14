import itertools
import re

from python_scripts.nanny.todo_utils import load_param
from python_scripts.nanny.xml_templates.hadrons_templates \
    import deep_copy_dict
from python_scripts.nanny.xml_templates.components import (
    gauge, eig, meson, highmodes, wrapper
)
import typing as t


def append_gauges_and_actions(modules: t.List, params: dict) -> None:

    modules += gauge.gauge_dp_params(
        ens=params['ENS'],
        series=params['SERIES'],
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

    for mass_label, mass in params['mass_dict'].items():

        modules += gauge.action_params(
            name=f"stag_{mass_label}",
            mass=mass,
            fat_label='gauge_fat',
            long_label='gauge_long',
            single_precision=False
        )

        modules += gauge.action_params(
            name=f"istag_{mass_label}",
            mass=mass,
            fat_label='gauge_fatf',
            long_label='gauge_longf',
            single_precision=True
        )


def append_epacks(modules: t.List, params: dict) -> None:

    modules += gauge.op_params('stag_op', 'stag_e')

    modules += eig.irl_params(
        name='epack',
        operator='stag_op_schur',
        lanczos_params={
            'alpha': params['ALPHA'],
            'beta': params['BETA'],
            'npoly': params['NPOLY'],
            'nstop': params['NSTOP'],
            'nk': params['NK'],
            'nm': params['NM']
        }
    )

    for mass_label, mass in params['mass_dict'].items():

        modules += eig.epack_modify_params(
            name=f"evecs_{mass_label}",
            epack="epack",
            mass=mass
        )


def append_mesons(modules: t.List,
                  operators: t.List[str], params: dict) -> None:

    mass_label = f"m{params['MASSES'][0]}"

    if 'local' in operators:
        modules += meson.meson_local_params(
            name=f"mf_ll_wv_local_{mass_label}",
            action=f"stag_{mass_label}",
            block=params['BLOCKSIZE'],
            lowmodes=f"evecs_{mass_label}",
            output=(f"e{params['EIGS']}n{params['NOISE']}dt{params['DT']}/"
                    f"mesons/{mass_label}/mf_{params['SERIES']}"),
            include_pion=True
        )

    if 'onelink' in operators:
        modules += meson.meson_onelink_params(
            name=f"mf_ll_wv_onelink_{mass_label}",
            action=f"stag_{mass_label}",
            block=params['BLOCKSIZE'],
            lowmodes=f"evecs_{mass_label}",
            gauge='gauge',
            output=(f"e{params['EIGS']}n{params['NOISE']}dt{params['DT']}/"
                    f"mesons/{mass_label}/mf_{params['SERIES']}")
        )


def append_solvers(modules: t.List, params: dict) -> None:

    for mass_label in params['mass_dict'].keys():

        modules += highmodes.mixed_precision_solver_params(
            name=f"stag_ama_{mass_label}",
            outer_action=f"stag_{mass_label}",
            inner_action=f"istag_{mass_label}",
            residual='1e-8'
        )

        modules += highmodes.ranLL_solver_params(
            name=f"stag_ranLL_{mass_label}",
            action=f"stag_{mass_label}",
            lowmodes=f"evecs_{mass_label}"
        )


def append_sources_and_quarks(modules: t.List, params: dict) -> None:

    for tslice in params['time_range']:
        modules += highmodes.source_params(
            name=f"noise_t{tslice}",
            nsrc=params['NOISE'],
            tstep=params['TIME'],
            t0=str(tslice)
        )

    def m1_eq_m2(x):
        return x[-2] == x[-1]

    for (tslice, glabel, slabel, mlabel, _) in filter(m1_eq_m2,
                                                      params['quark_iter']):

        quark = f"quark_{slabel}_{glabel}_{mlabel}_t{tslice}"
        source = f"noise_t{tslice}"
        solver = f"stag_{slabel}_{mlabel}"
        if slabel == "ama":
            guess = f"quark_ranLL_{glabel}_{mlabel}_t{tslice}"
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

    for (tslice, glabel, slabel, m1label, m2label) in filter(m1_ge_m2, params['quark_iter']):

        quark1 = f"quark_{slabel}_{glabel}_{m1label}_t{tslice}"
        quark2 = f"quark_{slabel}_pion_local_{m2label}_t{tslice}"

        if m1label == m2label:
            mass_label = m1label
        else:
            mass_label = m1label + "_" + m2label

        module_builder = getattr(highmodes, f"contractor_{glabel}_params")
        kwargs = {
            'name': f"corr_{slabel}_{glabel}_{mass_label}_t{tslice}",
            'source': quark1,
            'sink': quark2,
            'source_shift': f"noise_t{tslice}_shift",
            'output': (f"e{params['EIGS']}n{params['NOISE']}dt{params['DT']}/"
                       f"{mass_label}/{glabel}/{slabel}/corr_{glabel}"
                       f"_{slabel}_{mass_label}_t{tslice}_{params['SERIES']}")
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
        has_sea_mass = params['MASSES'][0] in x['name']
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
        lambda x: "_m" not in x['name'] and "_t" not in x['name']
    )

    sorted_modules = dp_gauges+indep_mass_tslice
    sorted_modules += meson_field_inputs+meson_fields
    sorted_modules += sp_gauges

    def gamma_order(x):
        for i, gamma in enumerate(['pion_local', 'vec_local', 'vec_onelink']):
            if gamma in x['name']:
                return i
        return -1

    def mass_order(x):
        for i, mass in enumerate(params['MASSES']):
            if mass in x['name']:
                return i
        return -1

    def mixed_mass_last(x):
        return len(re.findall(r'_m(\d+)', x['name']))

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
        runid=(f"LMI-RW-series-{params['SERIES']}"
               f"-{params['EIGS']}-eigs-{params['NOISE']}-noise"),
        sched=params.get('schedule', ''),
        cfg=params['CFG']
    )

    time = int(params["TIME"])
    t_start = int(params["TSTART"])
    t_stop = int(params["TSTOP"])+1
    t_step = int(params["DT"])
    params['time_range'] = list(range(t_start, t_stop, t_step))

    params["MASSES"] = params["MASSES"].strip().split(" ")
    mass_sufs = params["MASSES"]
    masses = [f"0.{suf}" for suf in mass_sufs]
    mass_labels = [f"m{suf}" for suf in mass_sufs]
    params['mass_dict'] = dict(zip(mass_labels, masses))

    gamma_labels = ['pion_local', 'vec_local', 'vec_onelink']
    solver_labels = ["ranLL", "ama"]
    params['quark_iter'] = list(itertools.product(
        params['time_range'],
        gamma_labels,
        solver_labels,
        mass_labels,
        mass_labels
    ))

    append_gauges_and_actions(modules, params)
    append_epacks(modules, params)

    meson_operators = tasks.get('meson', [])
    append_mesons(modules, meson_operators, params)

    if 'high_modes' in tasks:
        append_solvers(modules, params)
        append_sources_and_quarks(modules, params)
        append_contractors(modules, params)

    xml_params['grid']["modules"] = {"module": modules}

    module_info = [m["id"] for m in modules]
    schedule = build_schedule(module_info, params)

    return xml_params, schedule

    # return xml_params
