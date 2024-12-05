import copy
import os

def buildParams(**moduleTemplates):

    env = os.environ

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
                f"scheduleFile":"",
                "saveSchedule":"false",
                "parallelWriteMaxRetry":"-1",
            },
             "modules":{},
        },
    }

    modules = []

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge"
    module["options"]["file"] = f"lat/scidac/l{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge_fat"
    module["options"]["file"] = f"lat/scidac/fat{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge_long"
    module["options"]["file"] = f"lat/scidac/lng{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    params["grid"]["modules"] = {"module":modules}

    return params
