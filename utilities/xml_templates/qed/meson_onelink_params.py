import copy
import os

def build_params(**module_templates):

     params = {
         "grid":{
             "parameters":{
                 "runId":"LMI-RW-series-SERIES-EIGS-eigs-NOISE-noise",
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
                  f"schedule_file":"",
                 "saveSchedule":"false",
                 "parallelWriteMaxRetry":"-1",
             },
             "modules":{},
         },
     }

     modules = []

     module = copy.deepcopy(module_templates["load_gauge"])
     module["id"]["name"] = "gauge"
     module["options"]["file"] = "configs/lENSSERIES.ildg"
     modules.append(module)

     module = copy.deepcopy(module_templates["load_gauge"])
     module["id"]["name"] = "gauge_fat"
     module["options"]["file"] = "configs/fatENSSERIES.ildg"
     modules.append(module)

     module = copy.deepcopy(module_templates["load_gauge"])
     module["id"]["name"] = "gauge_long"
     module["options"]["file"] = "configs/lngENSSERIES.ildg"
     modules.append(module)

     module = copy.deepcopy(module_templates["epack_load"])
     module["id"]["name"] = "epack_l"
     module["options"]["filestem"] = "eigenpacks/eigENSnvSOURCEEIGSSERIES"
     module["options"]["size"] = "EIGS"
     modules.append(module)

     for i, m in enumerate(os.environ["MASSES"].strip().split(" ")):
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
          module["options"]["eigenPack"] = "epack_l"
          module["options"]["mass"] = mass
          modules.append(module)
          
          module = copy.deepcopy(module_templates["meson_field"])
          module["id"]["name"] = f"mf_ll_wv_{i}"
          module["options"].update({
               "action":f"stag_{i}",
               "block":"500",
               "spinTaste":{
                    "gammas":"(GX G1) (GY G1) (GZ G1)",
                    "gauge" :"gauge",
                    "applyG5":"false"
               },
               "lowModes":f"evecs_l_{i}",
               "output":f"eEIGSnNOISEdtDT/mesons/{mass_string}/mf_SERIES"
          })
          modules.append(module)

          break #only do one mass

     params["grid"]["modules"] = {"module":modules}

     return params
