from python_scripts.nanny.xml_templates.hadrons_templates import \
    module_templates
from python_scripts.utils import deep_copy_dict


def epack_load_params(name: str, filestem: str, eigs: str,
                      multifile: str = 'false') -> list:

    module = deep_copy_dict(module_templates["epack_load"])
    module["id"]["name"] = "epack"
    module["options"]["filestem"] = filestem
    module["options"]["size"] = str(eigs)
    module["options"]["multiFile"] = multifile

    return [module]


def eval_save_params(output: str) -> list:
    module = deep_copy_dict(module_templates["eval_save"])
    module["id"]["name"] = "eval_save"
    module["options"]["eigenPack"] = 'epack'
    module["options"]["output"] = output

    return [module]


def irl_params(name: str, operator: str, lanczos_params: dict,
               output: str = '', multifile: str = 'false',
               residual: str = '') -> list:

    module = deep_copy_dict(module_templates["irl"])
    module["id"]["name"] = "epack"
    module["options"]["op"] = "stag_op_schur"
    module["options"]["output"] = output
    module["options"]["multiFile"] = multifile

    lanczos_out = module["options"]["lanczosParams"]
    lanczos_out["Cheby"]["alpha"] = lanczos_params['alpha']
    lanczos_out["Cheby"]["beta"] = lanczos_params['beta']
    lanczos_out["Cheby"]["Npoly"] = lanczos_params['npoly']
    lanczos_out["Nstop"] = lanczos_params['nstop']
    lanczos_out["Nk"] = lanczos_params['nk']
    lanczos_out["Nm"] = lanczos_params['nm']

    try:
        _ = float(residual)
        lanczos_params["resid"] = residual
    except ValueError:
        pass

    return [module]


def epack_modify_params(name: str, epack: str, mass: str) -> list:
    module = deep_copy_dict(module_templates["epack_modify"])
    module["id"]["name"] = name
    module["options"]["eigenPack"] = epack
    module["options"]["mass"] = mass

    return [module]
