import copy
import os
import random

def build_params(**module_templates):

    jobid=int(random.random()*100)
    schedule_file=f"schedule/lma_meson_{jobid}.sched"
    runId=f"LMI-RW-series-{os.environ['SERIES']}-{os.environ['EIGS']}-eigs-{os.environ['NOISE']}-noise"
    random.seed(runId+os.environ["CFG"])
    gammas = {
        "pion"   :"(G5 G5)",
        "local"  :" ".join(["(GX GX)","(GY GY)","(GZ GZ)"]),
        "onelink": " ".join(["(GX G1)","(GY G1)","(GZ G1)"])
    }
    gammas_iter = list(gammas.items())

    # Make sure we iterate over pion first
    gammas_iter.sort(key=(lambda a: a[0] != "pion"))
    
    params = {
        "grid":{
            "parameters":{
                "runId":runId,
                "trajCounter":{
                    "start":"CFG",
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
                f"scheduleFile":schedule_file,
                "saveSchedule":"false",
                "parallelWriteMaxRetry":"-1",
            },
             "modules":{},
        },
    }

    modules = []

    module = copy.deepcopy(module_templates["load_gauge"])
    module["id"]["name"] = "gauge"
    module["options"]["file"] = "lat/scidac/lENSSERIES.ildg"
    modules.append(module)

    module = copy.deepcopy(module_templates["load_gauge"])
    module["id"]["name"] = "gauge_fat"
    module["options"]["file"] = "lat/scidac/fatENSSERIES.ildg"
    modules.append(module)

    module = copy.deepcopy(module_templates["load_gauge"])
    module["id"]["name"] = "gauge_long"
    module["options"]["file"] = "lat/scidac/lngENSSERIES.ildg"
    modules.append(module)

    module = copy.deepcopy(module_templates["cast_gauge"])
    module["id"]["name"] = "gauge_fatf"
    module["options"]["field"] = "gauge_fat"
    modules.append(module)

    module = copy.deepcopy(module_templates["cast_gauge"])
    module["id"]["name"] = "gauge_longf"
    module["options"]["field"] = "gauge_long"
    modules.append(module)

    module = copy.deepcopy(module_templates["sink"])
    module["id"]["name"] = "sink"
    module["options"]["mom"] = "0 0 0"
    modules.append(module)

    module = copy.deepcopy(module_templates["epack_load"])
    module["id"]["name"] = "epack"
    module["options"]["filestem"] = "eigen/eigENSnvSOURCEEIGSSERIES"
    module["options"]["size"] = "EIGS"
    modules.append(module)

    time = int(os.environ["TIME"])
    for block_index in range(int(os.environ["NBIAS"])):
        t0=random.randrange(time)
        block_label=f"n{block_index}"
        noise=f"noise_{block_label}"
        module = copy.deepcopy(module_templates["noise_rw"])
        module["id"]["name"] = noise
        module["options"]["nSrc"] = "1"
        module["options"]["tStep"] = str(time)
        module["options"]["t0"] = str(t0)
        #module["options"]["colorDiag"] = "false"
        modules.append(module)      

        for mass1_index, m1 in enumerate(os.environ["MASSES"].strip().split(" ")):
            mass1 = "0." + m1
            mass1_label = "m"+m1

            if block_index == 0:
                module = copy.deepcopy(module_templates["epack_modify"])
                module["id"]["name"] = f"evecs_{mass1_label}"
                module["options"]["eigenPack"] = "epack"
                module["options"]["mass"] = mass1
                modules.append(module)
                 
                module = copy.deepcopy(module_templates["action"])
                module["id"]["name"] = f"stag_{mass1_label}"
                module["options"]["mass"] = mass1
                module["options"]["gaugefat"] = "gauge_fat"
                module["options"]["gaugelong"] = "gauge_long"
                modules.append(module)

                module = copy.deepcopy(module_templates["action_float"])
                module["id"]["name"] = f"istag_{mass1_label}"
                module["options"]["mass"] = mass1
                module["options"]["gaugefat"] = "gauge_fatf"
                module["options"]["gaugelong"] = "gauge_longf"
                modules.append(module)
         
                module = copy.deepcopy(module_templates["mixed_precision_cg"])
                module["id"]["name"] = f"stag_fine_{mass1_label}"
                module["options"]["outerAction"] = f"stag_{mass1_label}"
                module["options"]["innerAction"] = f"istag_{mass1_label}"
                module["options"]["residual"] = "1e-8"
                modules.append(module)
          
                module = copy.deepcopy(module_templates["mixed_precision_cg"])
                module["id"]["name"] = f"stag_sloppy_{mass1_label}"
                module["options"]["outerAction"] = f"stag_{mass1_label}"
                module["options"]["innerAction"] = f"istag_{mass1_label}"
                module["options"]["residual"] = "2e-4"
                modules.append(module)
          
                module = copy.deepcopy(module_templates["lma_solver"])
                module["id"]["name"] = f"stag_ranLL_{mass1_label}"
                module["options"]["action"] = f"stag_{mass1_label}"
                module["options"]["lowModes"] = f"evecs_{mass1_label}"
                modules.append(module)

            for gamma_label, gamma_string in gammas_iter:

                # Only do sea mass for onelink
                if "onelink" in gamma_label and mass1_index != 0:
                    continue

                # Attach gauge link field to spin-taste for onelink
                gauge = "gauge" if gamma_label == "onelink" else ""

                for solver_label in ["ranLL","sloppy","fine"]:

                    solver=f"stag_{solver_label}_{mass1_label}"

                    if solver_label == "fine":
                        guess=f"quark_sloppy_{gamma_label}_{mass1_label}_{block_label}"
                    elif solver_label == "sloppy":
                        guess=f"quark_ranLL_{gamma_label}_{mass1_label}_{block_label}"
                    else:
                        guess=""
                        
                    quark_m1=f"quark_{solver_label}_{gamma_label}_{mass1_label}_{block_label}"
                    
                    module = copy.deepcopy(module_templates["quark_prop"])
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
                    for mass2_index, m2 in enumerate(os.environ["MASSES"].strip().split(" ")):

                        if mass2_index > mass1_index:
                            continue

                        mass2 = "0." + m2
                        mass2_label="m"+m2

                        if mass2_index == mass1_index:
                            mass_label = mass1_label
                        else:
                            mass_label = mass1_label + "_" + mass2_label

                        quark_m2=f"quark_{solver_label}_pion_{mass2_label}_{block_label}"

                        module = copy.deepcopy(module_templates["prop_contract"])
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
                            "output":f"eEIGSnNOISEdtDT/correlators/random_bias_colordiag/{mass_label}/{gamma_label}/bias/{solver_label}/corr_{gamma_label}_{solver_label}_{mass_label}_{block_label}_SERIES",
                        })
                        modules.append(module)

    params["grid"]["modules"] = {"module":modules}

    moduleList = [m["id"]["name"] for m in modules]

    f = open(schedule_file,"w")
    f.write(str(len(moduleList))+"\n"+"\n".join(moduleList))
    f.close()

    return params
