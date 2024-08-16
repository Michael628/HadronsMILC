import copy
import os
import re

def build_params(**module_templates):

    env = os.environ
    tStop = 1
    tStart = 0
    schedule_file=""

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
                "schedule_file":schedule_file,
                "saveSchedule":"false",
                "parallelWriteMaxRetry":"-1",
            },
            "modules":{},
        },
    }

    modules = []

    module = copy.deepcopy(module_templates["load_gauge"])
    module["id"]["name"] = "gauge_fat"
    module["options"]["file"] = f"configs/fat{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(module_templates["load_gauge"])
    module["id"]["name"] = "gauge"
    module["options"]["file"] = f"configs/l{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)
    
    module = copy.deepcopy(module_templates["load_gauge"])
    module["id"]["name"] = "gauge_long"
    module["options"]["file"] = f"configs/lng{env['ENS']}{env['SERIES']}.ildg"
    modules.append(module)

    module = copy.deepcopy(module_templates["epack_load"])
    module["id"]["name"] = "epack"
    module["options"]["filestem"] = f"eigenpacks/eig{env['ENS']}nv{env['SOURCEEIGS']}er8_grid_{env['SERIES']}"
    module["options"]["size"] = env['EIGS']
    module["options"]["multiFile"] = "true"
    modules.append(module)

    module = copy.deepcopy(module_templates["sink"])
    module["id"]["name"] = "sink"
    module["options"]["mom"] = "0 0 0"
    modules.append(module)

    time = int(env["TIME"])
    tStep = int(env["DT"])
    for time_index in range(tStart,tStop,tStep):
        block_label=f"t{time_index}"
        noise=f"noise_{block_label}"

        module = copy.deepcopy(module_templates["time_diluted_noise"])
        module["id"]["name"] = noise
        module["options"]["nsrc"] = env["NOISE"]
        module["options"]["tStep"] = str(time)
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
          
                module = copy.deepcopy(module_templates["lma_solver"])
                module["id"]["name"] = f"stag_ranLL_{mass1_label}"
                module["options"]["action"] = f"stag_{mass1_label}"
                module["options"]["lowModes"] = f"evecs_{mass1_label}"
                module["options"]["eigStart"] = '0'
                module["options"]["nEigs"] = env['EIGS']
                modules.append(module)

                module = copy.deepcopy(module_templates["meson_field"])
                module["id"]["name"] = f"mf_ll_wv_e_t0"
                module["options"].update({
                    "action":f"stag_{mass1_label}",
                    "block":env['BLOCKSIZE'],
                    "spinTaste":{
                        "gammas":"(G1 G1)",
                        "gauge" :"",
                        "applyG5":"false"
                    },
                    "right":noise+"_vec",
                    "lowModes":f"evecs_{mass1_label}",
                    "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/test/mesons/mf_e_{block_label}_{env['SERIES']}"
                })
                modules.append(module)

                module = copy.deepcopy(module_templates["meson_field"])
                module["id"]["name"] = f"mf_ll_wv_t0_e_local"
                module["options"].update({
                    "action":f"stag_{mass1_label}",
                    "block":env['BLOCKSIZE'],
                    "spinTaste":{
                        "gammas":"(G5 G5) (GX GX) (GY GY) (GZ GZ)",
                        "gauge" :"",
                        "applyG5":"false"
                    },
                    "left":noise+"_vec",
                    "lowModes":f"evecs_{mass1_label}",
                    "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/test/mesons/mf_{block_label}_e_{env['SERIES']}"
                })
                modules.append(module)

            for gamma_label, gamma_string in gammas_iter:

                # Only do sea mass for onelink
                if "onelink" in gamma_label and mass1_index != 0:
                    continue

                # Attach gauge link field to spin-taste for onelink
                gauge = "gauge" if "onelink" in gamma_label else ""

                for solver_label in ["ranLL"]:

                    solver=f"stag_{solver_label}_{mass1_label}"
                    guess=f"quark_ranLL_{gamma_label}_{mass1_label}_{block_label}" if "ama" in solver_label else ""
                    quark=f"quark_{solver_label}_{gamma_label}_{mass1_label}_{block_label}"

                    module = copy.deepcopy(module_templates["quark_prop"])
                    module["id"]["name"] = quark
                    module["options"].update({
                        "source"   :noise+"_vec",
                        "solver"   :solver,
                        "guess"    :guess,
                        "spinTaste":{
                            "gammas":gamma_string,
                            "gauge" :gauge,
                            "applyG5":"true"
                        }
                    })
                    modules.append(module)

                    quark_sink=f"quark_{solver_label}_pion_local_{mass1_label}_{block_label}"

                    module = copy.deepcopy(module_templates["prop_contract"])
                    module["id"]["name"] = f"corr_{solver_label}_{gamma_label}_{mass1_label}_{block_label}"
                    module["options"].update({
                        "source":quark,
                        "sink":quark_sink,
                        "sinkFunc":"sink",
                        "sourceShift":noise+"_shift",
                        "sourceGammas":gamma_string,
                        "sinkSpinTaste":{
                            "gammas":gamma_string,
                            "gauge" :gauge,
                            "applyG5":"true"
                        },
                        "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/test/{mass1_label}/{gamma_label}/{solver_label}/corr_{gamma_label}_{solver_label}_{mass1_label}_{block_label}_{env['SERIES']}",
                    })
                    modules.append(module)

    params["grid"]["modules"] = {"module":modules}

    moduleList = [m["id"]["name"] for m in modules]

    return params
