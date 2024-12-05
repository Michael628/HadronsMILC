import copy
import os
import re

def buildParams(**moduleTemplates):

    env = os.environ

    scheduleFile=f"schedules/lma2_{env['SERIES']}{env['CFG']}.sched"

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
                "scheduleFile":scheduleFile,
                "saveSchedule":"false",
                "parallelWriteMaxRetry":"-1",
            },
            "modules":{},
        },
    }

    modules = []

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge_fat"
    module["options"]["file"] = f"lat/scidac/fat{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge"
    module["options"]["file"] = f"lat/scidac/l{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)
    
    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge_long"
    module["options"]["file"] = f"lat/scidac/lng{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["epackLoad"])
    module["id"]["name"] = "epack"
    module["options"]["filestem"] = f"eigen/eig{env['ENS']}nv{env['SOURCEEIGS']}er8_grid_{env['SERIES']}"
    module["options"]["size"] = env['EIGS']
    module["options"]["multiFile"] = "true"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["sink"])
    module["id"]["name"] = "sink"
    module["options"]["mom"] = "0 0 0"
    modules.append(module)

    tStart = 0
    tStop = 1
    tStep = int(env["DT"])
    time = int(env['TIME'])
    for time_index in range(tStart,tStop,tStep):
        block_label=f"t{time_index}"
        noise=f"noise_{block_label}"

        module = copy.deepcopy(moduleTemplates["noiseRW"])
        module["id"]["name"] = noise
        module["options"]["nSrc"] = env["NOISE"]
        module["options"]["tStep"] = str(time)
        module["options"]["t0"] = str(time_index)
        modules.append(module)      

        for mass1_index, m1 in enumerate(masses):
            mass1 = "0." + m1
            mass1_label = "m"+m1

            if time_index == tStart:
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

                module = copy.deepcopy(moduleTemplates["castGauge"])
                module["id"]["name"] = "gauge_fatf"
                module["options"]["field"] = "gauge_fat"
                modules.append(module)

                module = copy.deepcopy(moduleTemplates["castGauge"])
                module["id"]["name"] = "gauge_longf"
                module["options"]["field"] = "gauge_long"
                modules.append(module)

                module = copy.deepcopy(moduleTemplates["actionF"])
                module["id"]["name"] = f"istag_{mass1_label}"
                module["options"]["mass"] = mass1
                module["options"]["gaugefat"] = "gauge_fatf"
                module["options"]["gaugelong"] = "gauge_longf"
                modules.append(module)
         
                module = copy.deepcopy(moduleTemplates["cgMP"])
                module["id"]["name"] = f"stag_ama_sloppy_{mass1_label}"
                module["options"]["outerAction"] = f"stag_{mass1_label}"
                module["options"]["innerAction"] = f"istag_{mass1_label}"
                module["options"]["residual"] = "1e-8"
                modules.append(module)
          
                module = copy.deepcopy(moduleTemplates["cgMP"])
                module["id"]["name"] = f"stag_ama_{mass1_label}"
                module["options"]["outerAction"] = f"stag_{mass1_label}"
                module["options"]["innerAction"] = f"istag_{mass1_label}"
                module["options"]["residual"] = "1e-10"
                modules.append(module)
          
                module = copy.deepcopy(moduleTemplates["lmaProj"])
                module["id"]["name"] = f"stag_ranLL_{mass1_label}"
                module["options"]["action"] = f"stag_{mass1_label}"
                module["options"]["lowModes"] = f"evecs_{mass1_label}"
                modules.append(module)

            for gamma_label, gamma_string in gammas_iter:

                # Only do sea mass for onelink
                if "onelink" in gamma_label and mass1_index != 0:
                    continue

                # Attach gauge link field to spin-taste for onelink
                gauge = "gauge" if "onelink" in gamma_label else ""

                for solver_label in ["ranLL","ama_sloppy",'ama']:

                    solver=f"stag_{solver_label}_{mass1_label}"
                    guess = ""
                    if solver_label == 'ama_sloppy':
                        guess=f"quark_ranLL_{gamma_label}_{mass1_label}_{block_label}"
                    if solver_label == "ama":
                        guess=f"quark_ama_sloppy_{gamma_label}_{mass1_label}_{block_label}"
                        
                    quark_m1=f"quark_{solver_label}_{gamma_label}_{mass1_label}_{block_label}"
                    
                    module = copy.deepcopy(moduleTemplates["quarkProp"])
                    module["id"]["name"] = quark_m1
                    module["options"].update({
                        "source"   :noise,
                        "solver"   :solver,
                        "guess"    :guess,
                        "spinTaste":{
                            "gammas":gamma_string,
                            "gauge" :gauge,
                            "applyG5":"true"
                        }
                    })
                    modules.append(module)

                    # Perform all cross-contractions between masses
                    for mass2_index, m2 in enumerate(masses):

                        if mass2_index > mass1_index:
                            continue

                        mass2 = "0." + m2
                        mass2_label="m"+m2

                        if mass2_index == mass1_index:
                            mass_label = mass1_label
                        else:
                            mass_label = mass1_label + "_" + mass2_label

                        quark_m2=f"quark_{solver_label}_pion_local_{mass2_label}_{block_label}"

                        module = copy.deepcopy(moduleTemplates["propContract"])
                        module["id"]["name"] = f"corr_{solver_label}_{gamma_label}_{mass_label}_{block_label}"
                        module["options"].update({
                            "source":quark_m1,
                            "sink":quark_m2,
                            "sinkFunc":"sink",
                            "sourceShift":noise+"_shift",
                            "sourceGammas":gamma_string,
                            "sinkSpinTaste":{
                                "gammas":gamma_string,
                                "gauge" :gauge,
                                "applyG5":"true"
                            },
                            "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/correlators/bias/{gamma_label}/{solver_label}/corr_{gamma_label}_{solver_label}_{mass_label}_{block_label}_{env['SERIES']}",
                        })
                        modules.append(module)

    params["grid"]["modules"] = {"module":modules}

    moduleList = [m["id"]["name"] for m in modules]

    f = open(scheduleFile,"w")
    f.write(str(len(moduleList))+"\n"+"\n".join(moduleList))
    f.close()
     
    return params
