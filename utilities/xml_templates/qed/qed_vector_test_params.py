import copy
import os
import random

def buildParams(**moduleTemplates):

    env = os.environ
    jobid=int(random.random()*100)
    scheduleFile=f"schedules/test_{env['SERIES']}{env['CFG']}.sched"
    gammas = {
        "pion_local"  :["(G5 G5)"],
        "vec_local"   :["(GX GX)","(GY GY)","(GZ GZ)"],
        #"vec_onelink" :["(GX G1)","(GY G1)","(GZ G1)"]
    }
    currents = {
        "current_local"  :" ".join(["(GX GX)","(GY GY)","(GZ GZ)","(GT GT)"]),
        #"current_onelink": " ".join(["(GX G1)","(GY G1)","(GZ G1)","(GT G1)"])
    }
    gammas_iter = list(gammas.items())

    params = {
        "grid":{
            "parameters":{
                "runId":f"A2A-series-{env['SERIES']}-{env['EIGS']}-eigs-{env['NOISE']}-noise",
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
                    
    module = copy.deepcopy(moduleTemplates["actionF"])
    module["id"]["name"] = f"istag_{mass_string}"
    module["options"]["mass"] = mass
    module["options"]["gaugefat"] = "gauge_fatf"
    module["options"]["gaugelong"] = "gauge_longf"
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

    seed=env["SEEDSTRING"]
    seed_string = f"{seed}0"

    module = copy.deepcopy(moduleTemplates["timeDilutedNoise"])
    module["id"]["name"] = f"ext_noise_{seed_string}"
    module["options"]["nsrc"] = env['NOISE']
    modules.append(module)
                    
    # Make external random wall sources
    for gamma_label, gamma_list in gammas_iter:
        for gamma in gamma_list:
            gamma_suf = gamma.replace("(","")
            gamma_suf = gamma_suf.replace(")","")
            gamma_suf = gamma_suf.replace(" ","_")
                            
            for solver_label in ["ranLL","ama"]:

                solver=f"stag_{solver_label}_{mass_string}"
                #solver+= "_subtract" if solver_label == "ama" else ""

                guess=f"quark_ranLL_{gamma_label}_{mass_string}_ext_{seed_string}_{gamma_suf}" if solver_label == "ama" else ""
                quark=f"quark_{solver_label}_{gamma_label}_{mass_string}_ext_{seed_string}_{gamma_suf}"
                    
                module = copy.deepcopy(moduleTemplates["quarkProp"])
                module["id"]["name"] = quark
                module["options"].update({
                    "source"   :f"ext_noise_{seed_string}_vec",
                    "solver"   :solver,
                    "guess"    :guess,
                    "spinTaste":{
                        "gammas":gamma,
                        "gauge" :"",
                        "applyG5":"false"
                    }
                })
                modules.append(module)

            module = copy.deepcopy(moduleTemplates["saveVector"])
            module["id"]["name"] = f"saveVecs_ext_{seed_string}_{gamma_suf}"
            module["options"]["field"] = f"quark_ama_{gamma_label}_{mass_string}_ext_{seed_string}_{gamma_suf}{gamma_suf}"
            module["options"]["output"] = f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/ext_vectors/{mass_string}/{gamma_suf}/{seed_string}_v"
            module["options"]["multiFile"] = "false"
            modules.append(module)

    for s in range(int(env["SEEDSTART"]),int(env["SEEDSTART"])+int(env["NSEEDS"])):
        seed_string=f"{seed}{s}"
        module = copy.deepcopy(moduleTemplates["timeDilutedNoise"])
        module["id"]["name"] = f"noise_{seed_string}"
        module["options"]["nsrc"] = env['NOISE']
        modules.append(module)
                    
        for solver_label in ["ranLL","ama"]:

            solver=f"stag_{solver_label}_{mass_string}"
            #solver+= "_subtract" if solver_label == "ama" else ""

            guess=f"quark_ranLL_{gamma_label}_{mass_string}_{seed_string}" if solver_label == "ama" else ""
            quark=f"quark_{solver_label}_{gamma_label}_{mass_string}_{seed_string}"
                    
            module = copy.deepcopy(moduleTemplates["quarkProp"])
            module["id"]["name"] = quark
            module["options"].update({
                "source"   :f"noise_{seed_string}_vec",
                "solver"   :solver,
                "guess"    :guess,
                "spinTaste":{
                    "gammas":"",
                    "gauge" :"",
                    "applyG5":"false"
                }
            })
            modules.append(module)

        module = copy.deepcopy(moduleTemplates["saveVector"])
        module["id"]["name"] = f"saveVecs_{seed_string}"
        module["options"]["field"] = f"quark_ama_{gamma_label}_{mass_string}_{seed_string}"
        module["options"]["output"] = f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/vectors/{mass_string}/{seed_string}_v"
        module["options"]["multiFile"] = "false"
        modules.append(module)

    params["grid"]["modules"] = {"module":modules}
     
    moduleList = [m["id"]["name"] for m in modules]

    f = open(scheduleFile, "w")
    f.write(str(len(moduleList)) + "\n" + "\n".join(moduleList))
    f.close()
    
    return params
