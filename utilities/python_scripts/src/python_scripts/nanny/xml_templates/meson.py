from python_scripts.nanny.xml_templates.hadrons_templates import \
    module_templates
from python_scripts.utils import deep_copy_dict


def meson_params(name: str, action: str, block: str,
                 lowmodes: str, gammas: str, gauge: str,
                 output: str) -> list:

    module = deep_copy_dict(module_templates["meson_field"])
    module["id"]["name"] = name
    module["options"]['action'] = action
    module["options"]['block'] = block
    module["options"]['lowModes'] = lowmodes
    module["options"]['spinTaste'] = {
        "gammas": gammas,
        "gauge": gauge,
        "applyG5": "false"
    }
    module["options"]['output'] = output

    return [module]


def meson_local_params(name: str, action: str, block: str,
                       lowmodes: str, output: str,
                       include_pion: bool = False) -> list:

    gammas = " ".join(["(GX GX)", "(GY GY)", "(GZ GZ)"])
    if include_pion:
        gammas = "(G5 G5) " + gammas

    return meson_params(name, action, block, lowmodes, gammas, '', output)


def meson_onelink_params(name: str, action: str, block: str,
                         lowmodes: str, gauge: str, output: str) -> list:
    assert gauge != ''

    gammas = " ".join(["(GX G1)", "(GY G1)", "(GZ G1)"])

    return meson_params(name, action, block, lowmodes, gammas, gauge, output)
