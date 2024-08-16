import copy
import os
import random

def buildParams(**moduleTemplates):

    env = os.environ

    scheduleFile=f"schedules/seq_sib_{env['SERIES']}{env['CFG']}.sched"

    masses=env["MASSES"].strip().split(" ")
    gammas = {
        "pion_local"  :["(G5 G5)"],
        "vec_local"   :["(GX GX)","(GY GY)","(GZ GZ)"],
        #"vec_onelink" :["(GX G1)","(GY G1)","(GZ G1)"]
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
                f"scheduleFile":scheduleFile,
                "saveSchedule":"false",
                "parallelWriteMaxRetry":"-1",
            },
             "modules":{},
        },
    }

    modules = []

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge"
    module["options"]["file"] = f"configs/l{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge_fat"
    module["options"]["file"] = f"configs/fat{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["loadGauge"])
    module["id"]["name"] = "gauge_long"
    module["options"]["file"] = f"configs/lng{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["castGauge"])
    module["id"]["name"] = "gauge_fatf"
    module["options"]["field"] = "gauge_fat"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["castGauge"])
    module["id"]["name"] = "gauge_longf"
    module["options"]["field"] = "gauge_long"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["sink"])
    module["id"]["name"] = "sink"
    module["options"]["mom"] = "0 0 0"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["epackLoad"])
    module["id"]["name"] = "epack"
    module["options"]["filestem"] = f"eigs/eig{env['ENS']}nv{env['SOURCEEIGS']}{env['SERIES']}"
    module["options"]["size"] = env['EIGS']
    module["options"]["multiFile"] = "false"
    modules.append(module)

    sib_label='sib'
    
    time = int(env["TIME"])
    tStart = int(env["TSTART"])
    tStop = int(env["TSTOP"])+1
    tStep = int(env["DT"])
    for time_index in range(tStart,tStop,tStep):
        block_label=f"t{time_index}"
        noise=f"noise_{block_label}"

        module = copy.deepcopy(moduleTemplates["noiseRW"])
        module["id"]["name"] = noise
        module["options"]["nSrc"] = env["NOISE"]
        module["options"]["tStep"] = str(time)
        module["options"]["t0"] = str(time_index)
        modules.append(module)      

        for mass_index, m1 in enumerate(masses):
            mass = "0." + m1
            mass_label = "m"+m1

            if time_index == tStart:
                module = copy.deepcopy(moduleTemplates["epackModify"])
                module["id"]["name"] = f"evecs_{mass_label}"
                module["options"]["eigenPack"] = "epack"
                module["options"]["mass"] = mass
                modules.append(module)
                 
                module = copy.deepcopy(moduleTemplates["action"])
                module["id"]["name"] = f"stag_{mass_label}"
                module["options"]["mass"] = mass
                module["options"]["gaugefat"] = "gauge_fat"
                module["options"]["gaugelong"] = "gauge_long"
                modules.append(module)

                module = copy.deepcopy(moduleTemplates["actionF"])
                module["id"]["name"] = f"istag_{mass_label}"
                module["options"]["mass"] = mass
                module["options"]["gaugefat"] = "gauge_fatf"
                module["options"]["gaugelong"] = "gauge_longf"
                modules.append(module)
         
                module = copy.deepcopy(moduleTemplates["cgMP"])
                module["id"]["name"] = f"stag_ama_{mass_label}"
                module["options"]["outerAction"] = f"stag_{mass_label}"
                module["options"]["innerAction"] = f"istag_{mass_label}"
                module["options"]["residual"] = "1e-8"
                modules.append(module)
          
                module = copy.deepcopy(moduleTemplates["lmaProj"])
                module["id"]["name"] = f"stag_ranLL_{mass_label}"
                module["options"]["action"] = f"stag_{mass_label}"
                module["options"]["lowModes"] = f"evecs_{mass_label}"
                module["options"]["eigStart"] = "0"
                module["options"]["nEigs"] = env['EIGS']
                modules.append(module)

            for gamma_label, gamma_list in gammas_iter:
                ext_gauge = "gauge" if "onelink" in gamma_label else ""

                for gamma in gamma_list:
                    gamma_suf = gamma.replace("(","")
                    gamma_suf = gamma_suf.replace(")","")
                    gamma_suf = gamma_suf.replace(" ","_")
                            
                    calc_pion = gamma_label == "pion_local"
                    
                    # Only do sea mass for onelink
                    if "onelink" in gamma_label and mass_index != 0:
                        continue

                    for solver_label in ["ranLL","ama"]:

                        solver=f"stag_{solver_label}_{mass_label}"
                        #guess=f"quark_ranLL_{gamma_label}_{mass_label}_{block_label}_{gamma_suf}" if solver_label == "ama" else ""
                        guess=""
                        quark_1=f"quark_{solver_label}_{gamma_label}_{mass_label}_{block_label}_{gamma_suf}"
                                
                        module = copy.deepcopy(moduleTemplates["quarkProp"])
                        module["id"]["name"] = quark_1
                        module["options"].update({
                            "source"   :noise,
                            "solver"   :solver,
                            "guess"    :guess,
                            "spinTaste":{
                                "gammas":gamma,
                                "gauge" :ext_gauge,
                                "applyG5": "false"
                            }
                        })
                        modules.append(module)

                        module = copy.deepcopy(moduleTemplates["seqGamma"])
                        module["id"]["name"] = quark_1 + f"_{sib_label}"
                        module["options"].update({
                            "q"   :quark_1+gamma_suf,
                            "tA"  :"0",
                            "tB"  :str(time),
                            "spinTaste":{
                                "gammas":"(G1 G1)",
                                "gauge" :"",
                                "applyG5": "false"
                            },
                            "mom":"0 0 0"
                        })
                        modules.append(module)
                            
                        #guess=f"quark_ranLL_{gamma_label}_{mass_label}_{block_label}_{gamma_suf}_{sib_label}_M" if solver_label == "ama" else ""
                        guess=""

                        module = copy.deepcopy(moduleTemplates["quarkProp"])
                        module["id"]["name"] = quark_1 + f"_{sib_label}_M"
                        module["options"].update({
                            "source"   :quark_1 + f"_{sib_label}",
                            "solver"   :solver,
                            "guess"    :guess,
                            "spinTaste":{
                                "gammas":"",
                                "gauge" :"",
                                "applyG5": "false"
                            }
                        })
                        modules.append(module)

                        if not calc_pion:

                            quark_1=quark_1+ f"_{sib_label}_M"
                            quark_2=f"quark_{solver_label}_pion_local_{mass_label}_{block_label}_G5_G5_{sib_label}_M"

                            module = copy.deepcopy(moduleTemplates["propContract"])
                            module["id"]["name"] = f"corr_{solver_label}_{gamma_label}_{mass_label}_{block_label}_{gamma_suf}_{sib_label}"
                            module["options"].update({
                                "source":quark_1,
                                "sink":quark_2,
                                "sinkFunc":"sink",
                                "sourceShift":noise+"_shift",
                                "sourceGammas":"",
                                "sinkSpinTaste":{
                                    "gammas":gamma,
                                    "gauge" :ext_gauge,
                                    "applyG5":"true"
                                },
                                "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/seq_correlators/{mass_label}/{gamma_label}_{sib_label}/{solver_label}/corr_{gamma_label}_{gamma_suf}_{sib_label}_{solver_label}_{mass_label}_{block_label}_{env['SERIES']}",
                            })
                            modules.append(module)

    params["grid"]["modules"] = {"module":modules}

    moduleList = [m["id"]["name"] for m in modules]

    f = open(scheduleFile,"w")
    f.write(str(len(moduleList))+"\n"+"\n".join(moduleList))
    f.close()
    
    return params
