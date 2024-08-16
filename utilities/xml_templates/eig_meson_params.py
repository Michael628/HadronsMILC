import copy
import os
import random

def build_params(**module_templates):

    env = os.environ
    jobid=int(random.random()*100)
    schedule_file=f"schedules/all_{jobid}.sched"

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

    module = copy.deepcopy(module_templates["action"])
    module["id"]["name"] = "stag_e"
    module["options"]["mass"] = "0.0"
    module["options"]["gaugefat"] = "gauge_fat"
    module["options"]["gaugelong"] = "gauge_long"
    modules.append(module)

    module = copy.deepcopy(module_templates["op"])
    module["id"]["name"] = "stag_op"
    module["options"]["action"] = "stag_e"
    modules.append(module)

    alpha = str(float( "%0.2g" % (float(env['ALPHA'])+float(env["BUNDLEINDEX"])*float(env["DALPHA"])))) if "DALPHA" in env.keys() else env['ALPHA']
    npoly = str(int(env['NPOLY'])+int(env["BUNDLEINDEX"])*int(env["DNPOLY"])) if "DNPOLY" in env.keys() else env['NPOLY']

    module = copy.deepcopy(module_templates["irl"])
    module["id"]["name"] = "epack"
    module["options"]["op"] = "stag_op_schur"
    module["options"]["lanczosParams"]["Cheby"]["alpha"] = alpha
    module["options"]["lanczosParams"]["Cheby"]["beta"] = env['BETA']
    module["options"]["lanczosParams"]["Cheby"]["Npoly"] = npoly
    module["options"]["lanczosParams"]["Nstop"] = env['NSTOP']
    module["options"]["lanczosParams"]["Nk"] = env['NK']
    module["options"]["lanczosParams"]["Nm"] = env['NM']
    if 'EIGOUT' in env.keys():
        module["options"]["output"] = env['EIGOUT'].format(ens=env['ENS'],series=env['SERIES'])
    if 'MULTIFILE' in env.keys():
        module["options"]["multiFile"] = env['MULTIFILE']
    if 'EIGRESID' in env.keys():
        module["options"]["lanczosParams"]["resid"] = env['EIGRESID']
    modules.append(module)
    
    module = copy.deepcopy(module_templates["sink"])
    module["id"]["name"] = "sink"
    module["options"]["mom"] = "0 0 0"
    modules.append(module)

    time = int(env["TIME"])
    tStart = int(env["TSTART"])
    tStop = int(env["TSTOP"])+1
    tStep = int(env["DT"])
    for time_index in range(tStart,tStop,tStep):
        block_label=f"t{time_index}"
        noise=f"noise_{block_label}"

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

                module = copy.deepcopy(module_templates["meson_field"])
                module["id"]["name"] = f"mf_ll_wv_local"
                module["options"].update({
                    "action":f"stag_{mass1_label}",
                    "block":"500",
                    "spinTaste":{
                        "gammas":gammas["pion_local"]+gammas["vec_local"],
                        "gauge" :"",
                        "applyG5":"false"
                    },
                    "lowModes":f"evecs_{mass1_label}",
                    "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons-dp/{mass1_label}/mf_{env['SERIES']}"
                })
                modules.append(module)

                module = copy.deepcopy(module_templates["meson_field"])
                module["id"]["name"] = f"mf_ll_wv_onelink"
                module["options"].update({
                    "action":f"stag_{mass1_label}",
                    "block":"250",
                    "spinTaste":{
                        "gammas":gammas["vec_onelink"],
                        "gauge" :"gauge",
                        "applyG5":"false"
                    },
                    "lowModes":f"evecs_{mass1_label}",
                    "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons-dp/{mass1_label}/mf_{env['SERIES']}"
                })
                modules.append(module)

    params["grid"]["modules"] = {"module":modules}

    moduleList = [m["id"]["name"] for m in modules]

    f = open(schedule_file,"w")
    f.write(str(len(moduleList))+"\n"+"\n".join(moduleList))
    f.close()
     
    return params
