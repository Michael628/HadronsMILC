import copy
import os
import random

def buildParams(**moduleTemplates):

    env = os.environ
    jobid=int(random.random()*100)
    scheduleFile=f"schedules/test_{jobid}.sched"

    params = {
        "grid":{
            "parameters":{
                "runId":f"{env['SEED']}-series-{env['SERIES']}-{env['EIGS']}-eigs-{env['NOISE']}-noise",
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
                "scheduleFile":scheduleFile,
                "saveSchedule":"false",
                "parallelWriteMaxRetry":"-1",
            },
            "modules":{},
        },
    }

    modules = []

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge"
    module["options"]["file"] = f'configs/l{env["ENS"]}{env["SERIES"]}.ildg'
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge_fat"
    module["options"]["file"] = f'configs/fat{env["ENS"]}{env["SERIES"]}.ildg'
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge_long"
    module["options"]["file"] = f'configs/lng{env["ENS"]}{env["SERIES"]}.ildg'
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["castGauge"])
    module["id"]["name"] = "gauge_fatf"
    module["options"]["field"] = "gauge_fat"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["castGauge"])
    module["id"]["name"] = "gauge_longf"
    module["options"]["field"] = "gauge_long"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["epackLoad"])
    module["id"]["name"] = "epack"
    module["options"]["filestem"] = f'eigs/eig{env["ENS"]}nv{env["SOURCEEIGS"]}{env["SERIES"]}'
    module["options"]["size"] = env["EIGS"]
    modules.append(module)

    m = os.environ["MASSES"].strip().split(" ")[0]
    mass_string = f"m{m}"
    mass = f"0.{m}"

    module = copy.deepcopy(moduleTemplates["action"])
    module["id"]["name"] = f"stag_{mass_string}"
    module["options"]["mass"] = mass
    module["options"]["gaugefat"] = "gauge_fat"
    module["options"]["gaugelong"] = "gauge_long"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["epackModify"])
    module["id"]["name"] = f"evecs_{mass_string}"
    module["options"]["eigenPack"] = "epack"
    module["options"]["mass"] = mass
    modules.append(module)
                    
    module = copy.deepcopy(moduleTemplates["timeDiluteNoise"])
    module["id"]["name"] = f"noise"
    module["options"]["nsrc"] = env['NOISE']
    modules.append(module)
                    
    module = copy.deepcopy(moduleTemplates["actionF"])
    module["id"]["name"] = f"istag_{mass_string}"
    module["options"]["mass"] = mass
    module["options"]["gaugefat"] = "gauge_fatf"
    module["options"]["gaugelong"] = "gauge_longf"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["sink"])
    module["id"]["name"] = "sink"
    module["options"]["mom"] = "0 0 0"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["cgMP"])
    module["id"]["name"] = f"stag_ama_{mass_string}"
    module["options"]["outerAction"] = f"stag_{mass_string}"
    module["options"]["innerAction"] = f"istag_{mass_string}"
    module["options"]["residual"] = "1e-8"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["lmaProj"])
    module["id"]["name"] = f"stag_ranLL_{mass_string}"
    module["options"]["action"] = f"stag_{mass_string}"
    module["options"]["lowModes"] = f"evecs_{mass_string}"
    modules.append(module)

    gamma_label="G1_G1"
    gamma="(G1 G1)"
    for solver_label in ["ranLL"]:#,"ama"]:

        solver=f"stag_{solver_label}_{mass_string}"
        solver+= "_subtract" if solver_label == "ama" else ""
        
        guess=f"quark_ranLL_{gamma_label}_{mass_string}" if solver_label == "ama" else ""
        quark=f"quark_{solver_label}_{gamma_label}_{mass_string}"
                    
        module = copy.deepcopy(moduleTemplates["quarkProp"])
        module["id"]["name"] = quark
        module["options"].update({
            "source"   :"noise_vec",
            "solver"   :solver,
            "guess"    :guess,
            "spinTaste":{
                "gammas":gamma,
                "gauge" :"",
                "applyG5":"false"
            }
        })
        modules.append(module)

    module = copy.deepcopy(moduleTemplates["a2aVector"])
    module["id"]["name"] = f"a2avec"
    module["options"]["noise"] = "noise"
    module["options"]["action"] = f"stag_{mass_string}"
    #module["options"]["lowModes"] = f"evecs_{mass_string}"
    module["options"]["lowModes"] = f"quark_ranLL_{gamma_label}_{mass_string}" + gamma_label
    module["options"]["solver"] = f"stag_ama_{mass_string}"
    module["options"]["highOutput"] = f"data/vectors/{env['SEED']}"
    module["options"]["highMultiFile"] = "true"
    module["options"]["norm2"] = env['NOISE']
    modules.append(module)

#        module = copy.deepcopy(moduleTemplates["propContract"])
#        module["id"]["name"] = f"corr_{solver_label}_{gamma_label}_{mass_string}"
#        module["options"].update({
#            "source":quark+gamma_label,
#            "sink":quark+gamma_label,
#            "sinkFunc":"sink",
#            "sourceShift":"noise_shift",
#            "sourceGammas":"",
#            "sinkSpinTaste":{
#                "gammas":gamma,
#                "gauge" :"",
#                "applyG5":"false"
#            },
#            "output":f"data/correlators/{mass_string}/{gamma_label}/{solver_label}/corr_{gamma_label}_{solver_label}_{mass_string}_{env['SERIES']}",
#        })
#        modules.append(module)
    
    quark=f"a2avec_v"
    solver_label="a2a"
    module = copy.deepcopy(moduleTemplates["propContract"])
    module["id"]["name"] = f"corr_{solver_label}_{gamma_label}_{mass_string}"
    module["options"].update({
        "source":quark,
        "sink":quark,
        "sinkFunc":"sink",
        "sourceShift":"noise_shift",
        "sourceGammas":"",
        "sinkSpinTaste":{
            "gammas":gamma,
            "gauge" :"",
            "applyG5":"false"
        },
        "output":f"data/correlators/{mass_string}/{gamma_label}/{solver_label}/corr_{gamma_label}_{solver_label}_{mass_string}_{env['SERIES']}",
    })
    modules.append(module)

    params["grid"]["modules"] = {"module":modules}
     
    moduleList = [m["id"]["name"] for m in modules]

    f = open(scheduleFile, "w")
    f.write(str(len(moduleList)) + "\n" + "\n".join(moduleList))
    f.close()

    return params
