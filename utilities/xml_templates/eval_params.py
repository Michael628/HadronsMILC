import copy
import os

def build_params(**module_templates):

    env = os.environ

    masses=env["MASSES"].strip().split(" ")
    gammas = {
        "pion_local"   :"(G5 G5)",
        "vec_local"  :" ".join(["(GX GX)","(GY GY)","(GZ GZ)"]),
        "vec_onelink": " ".join(["(GX G1)","(GY G1)","(GZ G1)"])
    }
    gammas_iter = list(gammas.items())

    # Make sure we iterate over pion first
    gammas_iter.sort(key=(lambda a: a[0] != "pion_local"))
    
    params = {
        "grid":{
            "parameters":{
                "runId":f"LMI-RW-series-{env['SERIES']}-{env['EIGS']}-eigs-{env['NOISE']}-noise",
                "trajCounter":{
                    "start":env["CFG"],
                    "end":"10000",
                    "step":"10000",
                },
                 "genetic":{
                     "popSize":"20",
                     "maxGen":"1000",
                     "maxCstGen":"100",
                     "mutationRate":"0.1",
                 },
                 "graphFile":"",
                f"schedule_file":"",
                "saveSchedule":"false",
                "parallelWriteMaxRetry":"-1",
            },
             "modules":{},
        },
    }

    modules = []

    module = copy.deepcopy(module_templates["epack_load"])
    module["id"]["name"] = "epack"
    module["options"]["filestem"] = f"eigs/eig{env['ENS']}nv{env['SOURCEEIGS']}{env['SERIES']}"
    module["options"]["size"] = env['EIGS']
    module["options"]["multiFile"] = "false"
    modules.append(module)

    module = copy.deepcopy(module_templates["eval_save"])
    module["id"]["name"] = "evals"
    module["options"]["eigenPack"] = "epack"
    module["options"]["output"] = f"eigs/evals/evalmassless{env['ENS']}nv{env['SOURCEEIGS']}{env['SERIES']}"
    modules.append(module)

    params["grid"]["modules"] = {"module":modules}

    moduleList = [m["id"]["name"] for m in modules]

    return params
