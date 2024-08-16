import copy
import os

def build_params(**module_templates):

     env = os.environ
     lat = env["ENS"][:4]
     
     gammas = ["(G5 G5)","(GX GX)","(GY GY)","(GZ GZ)"]
     gamma_string=" ".join(gammas)
     seeds = env["SEED"].split()

     params = {
         "grid":{
             "parameters":{
                 "runId":f"{env['SEED']}-series-{env['SERIES']}-{env['EIGS']}-eigs-{env['NOISE']}-noise",
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
                 "schedule_file":"",
                 "saveSchedule":"false",
                 "parallelWriteMaxRetry":"-1",
             },
             "modules":{},
         },
     }

     modules = []

     module = copy.deepcopy(module_templates["load_gauge"])
     module["id"]["name"] = "gauge"
     module["options"]["file"] = f'configs/l{env["ENS"]}{env["SERIES"]}.ildg'
     modules.append(module)

     module = copy.deepcopy(module_templates["load_gauge"])
     module["id"]["name"] = "gauge_fat"
     module["options"]["file"] = f'configs/fat{env["ENS"]}{env["SERIES"]}.ildg'
     modules.append(module)

     module = copy.deepcopy(module_templates["load_gauge"])
     module["id"]["name"] = "gauge_long"
     module["options"]["file"] = f'configs/lng{env["ENS"]}{env["SERIES"]}.ildg'
     modules.append(module)

     module = copy.deepcopy(module_templates["epack_load"])
     module["id"]["name"] = "epack_l"
     module["options"]["filestem"] = f'eigenpacks/eig{env["ENS"]}nv{env["SOURCEEIGS"]}{env["SERIES"]}'
     module["options"]["size"] = env["EIGS"]
     modules.append(module)

     m = os.environ["MASSES"].strip().split(" ")[0]
     mass_string = f"m{m}"
     mass = f"0.{m}"

     module = copy.deepcopy(module_templates["action"])
     module["id"]["name"] = f"stag_{mass_string}"
     module["options"]["mass"] = mass
     module["options"]["gaugefat"] = "gauge_fat"
     module["options"]["gaugelong"] = "gauge_long"
     modules.append(module)

     module = copy.deepcopy(module_templates["epack_modify"])
     module["id"]["name"] = f"evecs_l_{mass_string}"
     module["options"]["eigenPack"] = "epack_l"
     module["options"]["mass"] = mass
     modules.append(module)

     for i in range(len(seeds)):
         seedw = "w"+seeds[i]
         seedv = "v"+seeds[-1-i]

         module = copy.deepcopy(module_templates["load_vectors"])
         module["id"]["name"] = seedw
         module["options"]["filestem"] = f'data/vectors/{seedw[1:]}_w'
         #module["options"]["multiFile"] = 'true'
         module['options']['size'] = str(3*int(env['NOISE'])*int(env["TIME"]))
         modules.append(module)
     
         module = copy.deepcopy(module_templates["load_vectors"])
         module["id"]["name"] = seedv
         module["options"]["filestem"] = f'data/vectors/{seedv[1:]}_v'
         #module["options"]["multiFile"] = 'true'
         module['options']['size'] = str(3*int(env['NOISE'])*int(env["TIME"]))
         modules.append(module)
     
         module = copy.deepcopy(module_templates["meson_field"])
         module["id"]["type"] = "MContraction::QEDMesonField"
         module["id"]["name"] = f"mf_ll_{seedw}_eig"
         module["options"].update({
             "action":f"stag_{mass_string}",
             "block":"300",
             "left":seedw,
             "right":"",
             "em_func":"",
             "EmSeedString":"test_1",
             "nem_fields":"0",
             "spinTaste":{
                 "gammas":gamma_string,
                 "gauge" :"gauge" if "G1" in gamma_string else "",
                 "applyG5":"false"
             },
             "lowModes":f"evecs_l_{mass_string}",
             "output":f"data/mesons/mf_{env['SERIES']}_{seedw}_eig{env['EIGS']}"
         })
         modules.append(module)                    

         module = copy.deepcopy(module_templates["meson_field"])
         module["id"]["type"] = "MContraction::QEDMesonField"
         module["id"]["name"] = f"mf_ll_eig_{seedv}"
         module["options"].update({
             "action":f"stag_{mass_string}",
             "block":"300",
             "left":"",
             "right":seedv,
             "em_func":"",
             "EmSeedString":"test_1",
             "nem_fields":"0",
             "spinTaste":{
                 "gammas":gamma_string,
                 "gauge" :"gauge" if "G1" in gamma_string else "",
                 "applyG5":"false"
             },
             "lowModes":f"evecs_l_{mass_string}",
             "output":f"data/mesons/mf_{env['SERIES']}_eig{env['EIGS']}_{seedv}"
         })
         modules.append(module)                    

         module = copy.deepcopy(module_templates["meson_field"])
         module["id"]["type"] = "MContraction::QEDMesonField"
         module["id"]["name"] = f"mf_ll_{seedw}_{seedv}"
         module["options"].update({
              "action":f"stag_{mass_string}",
              "block":"300",
              "left":seedw,
              "right":seedv,
              "em_func":"",
              "EmSeedString":"test_1",
              "nem_fields":"0",
              "spinTaste":{
                   "gammas":gamma_string,
                   "gauge" :"gauge" if "G1" in gamma_string else "",
                   "applyG5":"false"
              },
              "lowModes":f"",
              "output":f"data/mesons/mf_{env['SERIES']}_{seedw}_{seedv}"
         })
         modules.append(module)                    
         
     params["grid"]["modules"] = {"module":modules}
     
     moduleList = [m["id"]["name"] for m in modules]

     return params
