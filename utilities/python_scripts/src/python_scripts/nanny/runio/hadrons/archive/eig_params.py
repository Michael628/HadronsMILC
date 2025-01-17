import copy
import os


def build_params(**module_templates):

    env = os.environ

    params = {
        "grid": {
            "parameters": {
                "runId": f"LMI-RW-series-{env['SERIES']}-{env['EIGS']}-eigs-{env['NOISE']}-noise",
                "trajCounter": {
                    "start": env['CFG'],
                    "end": "10000",
                    "step": "10000",
                },
                "genetic": {
                    "popSize": "20",
                    "maxGen": "1000",
                    "maxCstGen": "100",
                    "mutationRate": "0.1",
                },
                "graphFile": "",
                "scheduleFile": "",
                "saveSchedule": "false",
                "parallelWriteMaxRetry": "-1",
            },
            "modules": {},
        },
    }

    modules = []

    module = copy.deepcopy(module_templates["load_gauge"])
    module["id"]["name"] = "gauge_fat"
    module["options"]["file"] = f"lat/scidac/fat{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(module_templates["load_gauge"])
    module["id"]["name"] = "gauge_long"
    module["options"]["file"] = f"lat/scidac/lng{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(module_templates["action"])
    module["id"]["name"] = "stag_e"
    module["options"]["mass"] = "0.0"
    module["options"]["gaugefat"] = "gauge_fat"
    module["options"]["gaugelong"] = "gauge_long"
    modules.append(module)

    module = copy.deepcopy(module_templates["op"])
    module["id"]["name"] = "stag_op"
    module["options"]["action"] = "stag_e"
    modules.append(module)

    alpha = str(float("%0.2g" % (float(env['ALPHA'])+float(env["BUNDLEINDEX"])*float(
        env["DALPHA"])))) if "DALPHA" in env.keys() else env['ALPHA']
    npoly = str(int(env['NPOLY'])+int(env["BUNDLEINDEX"]) *
                int(env["DNPOLY"])) if "DNPOLY" in env.keys() else env['NPOLY']

    module = copy.deepcopy(module_templates["irl"])
    module["id"]["name"] = "epack"
    module["options"]["op"] = "stag_op_schur"
    module["options"]["lanczosParams"]["Cheby"]["alpha"] = alpha
    module["options"]["lanczosParams"]["Cheby"]["beta"] = env['BETA']
    module["options"]["lanczosParams"]["Cheby"]["Npoly"] = npoly
    module["options"]["lanczosParams"]["Nstop"] = env['NSTOP']
    module["options"]["lanczosParams"]["Nk"] = env['NK']
    module["options"]["lanczosParams"]["Nm"] = env['NM']
    if 'EIGOUT' in env.keys():
        module["options"]["output"] = env['EIGOUT'].format(
            ens=env['ENS'], series=env['SERIES'])
    if 'MULTIFILE' in env.keys():
        module["options"]["multiFile"] = env['MULTIFILE']
    if 'EIGRESID' in env.keys():
        module["options"]["lanczosParams"]["resid"] = env['EIGRESID']
    modules.append(module)

    params["grid"]["modules"] = {"module": modules}

    return params
