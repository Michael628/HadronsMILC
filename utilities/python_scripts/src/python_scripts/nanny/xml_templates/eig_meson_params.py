import copy
import os
import random

def buildParams(**moduleTemplates):

    env = os.environ
    jobid=int(random.random()*100)
    scheduleFile=f"schedules/eig_meson_{jobid}.sched"

    masses=env["MASSES"].strip().split(" ")

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

    module = copy.deepcopy(moduleTemplates["action"])
    module["id"]["name"] = "stag_e"
    module["options"]["mass"] = "0.0"
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

    module = copy.deepcopy(moduleTemplates["op"])
    module["id"]["name"] = "stag_op"
    module["options"]["action"] = "stag_e"
    modules.append(module)

    module = copy.deepcopy(moduleTemplates["IRL"])
    module["id"]["name"] = "epack"
    module["options"]["op"] = "stag_op_schur"
    module["options"]["lanczosParams"]["Cheby"]["alpha"] = env['ALPHA']
    module["options"]["lanczosParams"]["Cheby"]["beta"] = env['BETA']
    module["options"]["lanczosParams"]["Cheby"]["Npoly"] = env['NPOLY']
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
    
    module = copy.deepcopy(moduleTemplates["sink"])
    module["id"]["name"] = "sink"
    module["options"]["mom"] = "0 0 0"
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

        if mass1_index == 0:
            module = copy.deepcopy(moduleTemplates["mesonField"])
            module["id"]["name"] = f"mf_ll_wv_onelink"
            module["options"].update({
                "action":f"stag_{mass1_label}",
                "block":env['BLOCKSIZE'],
                "spinTaste":{
                    "gammas":"(GX G1) (GY G1) (GZ G1)",
                    "gauge" :"gauge",
                    "applyG5":"false"
                },
                "lowModes":f"evecs_{mass1_label}",
                "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons/{(mass1_label+'/') if len(masses) > 1 else ''}mf_{env['SERIES']}"
            })
            modules.append(module)

            module = copy.deepcopy(moduleTemplates["mesonField"])
            module["id"]["name"] = f"mf_ll_wv_local"
            module["options"].update({
                "action":f"stag_{mass1_label}",
                "block":env['BLOCKSIZE'],
                "spinTaste":{
                    "gammas":"(G5 G5) (GX GX) (GY GY) (GZ GZ)",
                    "gauge" :"",
                    "applyG5":"false"
                },
                "lowModes":f"evecs_{mass1_label}",
                "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons/{(mass1_label+'/') if len(masses) > 1 else ''}mf_{env['SERIES']}"
            })
            modules.append(module)

    params["grid"]["modules"] = {"module":modules}

    moduleList = [m["id"]["name"] for m in modules]

    f = open(scheduleFile,"w")
    f.write(str(len(moduleList))+"\n"+"\n".join(moduleList))
    f.close()
     
    return params
