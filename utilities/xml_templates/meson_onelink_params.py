import copy
import os

def buildParams(**moduleTemplates):

    env = os.environ

    masses=env["MASSES"].strip().split(" ")

    gammas = {
         "pion_local"   :"(G5 G5)",
         "vec_local"  :" ".join(["(GX GX)","(GY GY)","(GZ GZ)"]),
         "vec_onelink": " ".join(["(GX G1)","(GY G1)","(GZ G1)"])
    }
    gammas_iter = list(gammas.items())

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
                "scheduleFile":"",
                "saveSchedule":"false",
                "parallelWriteMaxRetry":"-1",
            },
            "modules":{},
        },
    }

    modules = []

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge_fat"
    module["options"]["file"] = f"configs/fat{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge"
    module["options"]["file"] = f"configs/l{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)
    
    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge_long"
    module["options"]["file"] = f"configs/lng{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["action"])
    module["id"]["name"] = "stag_e"
    module["options"]["mass"] = "0.0"
    module["options"]["gaugefat"] = "gauge_fat"
    module["options"]["gaugelong"] = "gauge_long"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["epackLoad"])
    module["id"]["name"] = "epack"
    module["options"]["filestem"] = f"eigs/eig{env['ENS']}nv{env['SOURCEEIGS']}{env['SERIES']}"
    module["options"]["size"] = env['EIGS']
    module["options"]["multiFile"] = "false"
    modules.append(module)

    for mass1_index, m1 in enumerate(masses):
         mass1 = "0." + m1
         mass1_label = "m"+m1

         module = copy.deepcopy(moduleTemplates["epackModify"])
         module["id"]["name"] = f"evecs_{mass1_label}"
         module["options"]["eigenPack"] = "epack"
         module["options"]["mass"] = mass1
         modules.append(module)
                 
         module = copy.deepcopy(moduleTemplates["action"])
         module["id"]["name"] = f"stag_{mass1_label}"
         module["options"]["mass"] = mass1
         module["options"]["gaugefat"] = "gauge_fat"
         module["options"]["gaugelong"] = "gauge_long"
         modules.append(module)
         
         module = copy.deepcopy(moduleTemplates["mesonField"])
         module["id"]["name"] = f"mf_ll_wv_onelink"
         module["options"].update({
              "action":f"stag_{mass1_label}",
              "block":"500",
              "spinTaste":{
                   "gammas":gammas["vec_onelink"],
                   "gauge" :"gauge",
                   "applyG5":"false"
              },
              "lowModes":f"evecs_{mass1_label}",
              "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons/{mass1_label}/mf_{env['SERIES']}"
         })
         modules.append(module)
         
         break
           
    params["grid"]["modules"] = {"module":modules}

    return params
