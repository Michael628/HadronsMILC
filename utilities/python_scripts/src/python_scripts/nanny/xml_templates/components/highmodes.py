from python_scripts.nanny.xml_templates.hadrons_templates import \
    module_templates
from python_scripts.utils import deep_copy_dict


# Sources
def source_params(name: str, nsrc: str, tstep: str, t0: str) -> list:

    module = deep_copy_dict(module_templates["noise_rw"])
    module["id"]["name"] = name
    module["options"]["nSrc"] = nsrc
    module["options"]["tStep"] = tstep
    module["options"]["t0"] = t0

    return [module]


# Solvers
def mixed_precision_solver_params(name: str, outer_action: str,
                                  inner_action: str, residual: str) -> list:

    module = deep_copy_dict(module_templates["mixed_precision_cg"])
    module["id"]["name"] = name
    module["options"]["outerAction"] = outer_action
    module["options"]["innerAction"] = inner_action
    module["options"]["residual"] = residual

    return [module]


def ranLL_solver_params(name: str, action: str, lowmodes: str):

    module = deep_copy_dict(module_templates["lma_solver"])
    module["id"]["name"] = name
    module["options"]["action"] = action
    module["options"]["lowModes"] = lowmodes

    return [module]


# Propagators
def propagator_params(name: str, source: str, solver: str,
                      guess: str, gammas: str, gauge: str):
    module = deep_copy_dict(module_templates["quark_prop"])
    module["id"]["name"] = name
    module["options"].update({
        "source": source,
        "solver": solver,
        "guess": guess,
        "spinTaste": {
            "gammas": gammas,
            "gauge": gauge,
            "applyG5": "true"
        }
    })

    return [module]


def propagator_pion_local_params(name: str, source: str, solver: str,
                                 guess: str):
    return propagator_params(
        name=name,
        source=source,
        solver=solver,
        guess=guess,
        gammas="(G5 G5)",
        gauge=''
    )


def propagator_vec_local_params(name: str, source: str, solver: str,
                                guess: str):
    return propagator_params(
        name=name,
        source=source,
        solver=solver,
        guess=guess,
        gammas="(GX GX) (GY GY) (GZ GZ)",
        gauge=''
    )


def propagator_vec_onelink_params(name: str, source: str, solver: str,
                                  guess: str, gauge: str):
    return propagator_params(
        name=name,
        source=source,
        solver=solver,
        guess=guess,
        gammas="(GX G1) (GY G1) (GZ G1)",
        gauge=gauge
    )


# Contractors
def sink_params():
    module = deep_copy_dict(module_templates["sink"])
    module["id"]["name"] = "sink"
    module["options"]["mom"] = "0 0 0"
    return [module]


def contractor_params(name: str, source: str, source_shift: str, sink: str,
                      source_gammas: str, sink_gammas: str,
                      sink_gauge: str, output: str):
    module = deep_copy_dict(module_templates["prop_contract"])
    module["id"]["name"] = name
    module["options"].update({
        "source": source,
        "sourceShift": source_shift,
        "sink": sink,
        "sinkFunc": "sink",
        "sourceGammas": source_gammas,
        "sinkSpinTaste": {
            "gammas": sink_gammas,
            "gauge": sink_gauge,
            "applyG5": "true"
        },
        "output": output
    })

    return [module]


def contractor_pion_local_params(name: str, source: str, source_shift: str,
                                 sink: str, output: str):
    return contractor_params(
        name=name,
        source=source,
        source_shift=source_shift,
        sink=sink,
        source_gammas="(G5 G5)",
        sink_gammas="(G5 G5)",
        sink_gauge='',
        output=output
    )


def contractor_vec_local_params(name: str, source: str, source_shift: str,
                                sink: str, output: str):
    return contractor_params(
        name=name,
        source=source,
        source_shift=source_shift,
        sink=sink,
        source_gammas="(GX GX) (GY GY) (GZ GZ)",
        sink_gammas="(GX GX) (GY GY) (GZ GZ)",
        sink_gauge='',
        output=output
    )


def contractor_vec_onelink_params(name: str, source: str, source_shift: str,
                                  sink: str, gauge: str, output: str):
    return contractor_params(
        name=name,
        source=source,
        source_shift=source_shift,
        sink=sink,
        source_gammas="(GX G1) (GY G1) (GZ G1)",
        sink_gammas="(GX G1) (GY G1) (GZ G1)",
        sink_gauge=gauge,
        output=output
    )
