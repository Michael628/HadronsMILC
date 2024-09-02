import copy
import os

def build_params(**module_templates):

     env = os.environ
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

     module = copy.deepcopy(module_templates["epack_load"])
     module["id"]["name"] = "epack"
     module["options"]["filestem"] = f"eigen/eig{env['ENS']}nv{env['SOURCEEIGS']}er8_grid_{env['SERIES']}"
     module["options"]["size"] = env['EIGS']
     module["options"]["multiFile"] = "true"
     modules.append(module)

     for i, m in enumerate(masses):
          mass = "0." + m
          mass_string = "m"+mass[2:]
          module = copy.deepcopy(module_templates["action"])
          module["id"]["name"] = f"stag_{i}"
          module["options"]["mass"] = mass
          module["options"]["gaugefat"] = "gauge_fat"
          module["options"]["gaugelong"] = "gauge_long"
          modules.append(module)

          module = copy.deepcopy(module_templates["epack_modify"])
          module["id"]["name"] = f"evecs_l_{i}"
          module["options"]["eigenPack"] = "epack"
          module["options"]["mass"] = mass
          modules.append(module)
          
          module = copy.deepcopy(module_templates["meson_field"])
          module["id"]["name"] = f"mf_ll_wv_onelink_{i}"
          module["options"].update({
               "action":f"stag_{i}",
               "block":"200",
               "spinTaste":{
                    "gammas":"(GX G1) (GY G1) (GZ G1)",
                    "gauge" :"gauge",
                    "applyG5":"false"
               },
               "lowModes":f"evecs_l_{i}",
               "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons/{(mass_string+'/') if len(masses) > 1 else ''}mf_{env['SERIES']}"
          })
          modules.append(module)

          module = copy.deepcopy(module_templates["meson_field"])
          module["id"]["name"] = f"mf_ll_wv_local_{i}"
          module["options"].update({
               "action":f"stag_{i}",
               "block":"200",
               "spinTaste":{
                    "gammas":"(G5 G5) (GX GX) (GY GY) (GZ GZ)",
                    "gauge" :"",
                    "applyG5":"false"
               },
               "lowModes":f"evecs_l_{i}",
               "output":f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/mesons/{(mass_string+'/') if len(masses) > 1 else ''}mf_{env['SERIES']}"
          })
          modules.append(module)

          break #only do one mass

     params["grid"]["modules"] = {"module":modules}

     return params
