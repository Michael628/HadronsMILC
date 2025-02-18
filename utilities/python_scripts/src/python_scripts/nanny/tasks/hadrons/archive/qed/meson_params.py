import copy
import os

def build_params(**module_templates):

    env = os.environ

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
                f"scheduleFile":"",
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

    for mass1_index, m1 in enumerate(masses):
         mass1 = "0." + m1
         mass1_label = "m"+m1

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

    params["grid"]["modules"] = {"module":modules}

    moduleList = [m["id"]["name"] for m in modules]
    
    return params
