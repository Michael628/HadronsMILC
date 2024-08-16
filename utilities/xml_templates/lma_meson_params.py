import copy
import os
import random

def build_params(**module_templates):

    env = os.environ

    jobid=int(random.random()*100)
    schedule_file=f"schedules/lma_{jobid}.sched"
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
                f"schedule_file":schedule_file,
                "saveSchedule":"false",
                "parallelWriteMaxRetry":"-1",
            },
             "modules":{},
        },
    }

    modules = []

    module = copy.deepcopy(module_templates["load_gauge"])
    module["id"]["name"] = "gauge"
    module["options"]["file"] = f"lat/scidac/l{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(module_templates["load_gauge"])
    module["id"]["name"] = "gauge_fat"
    module["options"]["file"] = f"lat/scidac/fat{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(module_templates["load_gauge"])
    module["id"]["name"] = "gauge_long"
    module["options"]["file"] = f"lat/scidac/lng{env['ENS']}{env['SERIES']}.ildg"
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
    module["options"]["filestem"] = f"eigen/eig{env['ENS']}nv{env['SOURCEEIGS']}er8_grid_{env['SERIES']}"
    module["options"]["size"] = env['EIGS']
    module["options"]["multiFile"] = "true"
    modules.append(module)

    time = int(env["TIME"])
    tStart = int(env["TSTART"])
    tStop = int(env["TSTOP"])+1
    tStep = int(env["DT"])
    for time_index in range(tStart,tStop,tStep):
        block_label=f"t{time_index}"
        noise=f"noise_{block_label}"

        module = copy.deepcopy(module_templates["noise_rw"])
        module["id"]["name"] = noise
        module["options"]["nSrc"] = env["NOISE"]
        module["options"]["tStep"] = str(time)
        module["options"]["t0"] = str(time_index)
        modules.append(module)      

        for mass1_index, m1 in enumerate(masses):
            mass1 = "0." + m1
            mass1_label = "m"+m1

            if time_index == tStart:
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

                if mass1_index == 0:
                    module = copy.deepcopy(module_templates["meson_field"])
                    module["id"]["name"] = f"mf_ll_wv_onelink"
                    module["options"].update({
                        "action":f"stag_{mass1_label}",
                        "block":"200",
                        "spinTaste":{
                            "gammas":"(GX G1) (GY G1) (GZ G1)",
                            "gauge" :"gauge",
                            "applyG5":"false"
                        },
                        "lowModes":f"evecs_{mass1_label}",
                        "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons/{(mass1_label+'/') if len(masses) > 1 else ''}mf_{env['SERIES']}"
                    })
                    modules.append(module)

                    module = copy.deepcopy(module_templates["meson_field"])
                    module["id"]["name"] = f"mf_ll_wv_local"
                    module["options"].update({
                        "action":f"stag_{mass1_label}",
                        "block":"200",
                        "spinTaste":{
                            "gammas":"(G5 G5) (GX GX) (GY GY) (GZ GZ)",
                            "gauge" :"",
                            "applyG5":"false"
                        },
                        "lowModes":f"evecs_{mass1_label}",
                        "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons/{(mass1_label+'/') if len(masses) > 1 else ''}mf_{env['SERIES']}"
                    })
                    modules.append(module)

                module = copy.deepcopy(module_templates["action_float"])
                module["id"]["name"] = f"istag_{mass1_label}"
                module["options"]["mass"] = mass1
                module["options"]["gaugefat"] = "gauge_fatf"
                module["options"]["gaugelong"] = "gauge_longf"
                modules.append(module)
         
                module = copy.deepcopy(module_templates["mixed_precision_cg"])
                module["id"]["name"] = f"stag_ama_{mass1_label}"
                module["options"]["outerAction"] = f"stag_{mass1_label}"
                module["options"]["innerAction"] = f"istag_{mass1_label}"
                module["options"]["residual"] = "1e-8"
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
                gauge = "gauge" if "onelink" in gamma_label else ""

                for solver_label in ["ranLL","ama"]:

                    solver=f"stag_{solver_label}_{mass1_label}"
                    guess=f"quark_ranLL_{gamma_label}_{mass1_label}_{block_label}" if solver_label == "ama" else ""
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
                            "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/correlators/{mass_label+'/' if len(masses) > 1 else ''}{gamma_label}/{solver_label}/corr_{gamma_label}_{solver_label}_{mass_label}_{block_label}_{env['SERIES']}",
                        })
                        modules.append(module)

    params["grid"]["modules"] = {"module":modules}

    moduleList = [m["id"]["name"] for m in modules]
    moduleList.sort(key=lambda x: 1 if "mf" in x else 0)
    
    f = open(schedule_file,"w")
    f.write(str(len(moduleList))+"\n"+"\n".join(moduleList))
    f.close()
    
    return params
