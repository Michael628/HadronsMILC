import copy
import os
def build_params(**module_templates):

     env = os.environ

     params = {
         "grid":{
             "parameters":{
                 "runId":f"LMI-RW-series-{env['SERIES']}-{env['EIGS']}-eigs-{env['NOISE']}-noise",
                 "trajCounter":{
                     "start":env['CFG'],
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
     module["id"]["name"] = "gauge_fat"
     module["options"]["file"] = f"configs/fat{env['ENS']}{env['SERIES']}.ildg"
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

     module = copy.deepcopy(module_templates["irl"])
     module["id"]["name"] = "epack_l"
     module["options"]["op"] = "stag_op_schur"
     module["options"]["lanczosParams"]["Cheby"]["alpha"] = "0.01"
     module["options"]["lanczosParams"]["Cheby"]["beta"] = "24"
     module["options"]["lanczosParams"]["Cheby"]["Npoly"] = "101"
     module["options"]["lanczosParams"]["Nstop"] = "1000"
     module["options"]["lanczosParams"]["Nk"] = "1050"
     module["options"]["lanczosParams"]["Nm"] = "1300"
     module["options"]["output"] = f"eigs/local/eig{env['ENS']}nv2000er8_grid_{env['SERIES']}"
     module["options"]["multiFile"] = "false"
     modules.append(module)

     params["grid"]["modules"] = {"module":modules}

     return params
