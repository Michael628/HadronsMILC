import copy
import os
def buildParams(**moduleTemplates):

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
                 "scheduleFile":"",
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
     module["id"]["name"] = "gauge_long"
     module["options"]["file"] = f"lat/scidac/lng{env['ENS']}{env['SERIES']}.ildg"
     modules.append(module)

     module = copy.deepcopy(moduleTemplates["action"])
     module["id"]["name"] = "stag_e"
     module["options"]["mass"] = "0.0"
     module["options"]["gaugefat"] = "gauge_fat"
     module["options"]["gaugelong"] = "gauge_long"
     modules.append(module)

     
     module = copy.deepcopy(moduleTemplates["op"])
     module["id"]["name"] = "stag_op"
     module["options"]["action"] = "stag_e"
     modules.append(module)

     alpha = str(float( "%0.2g" % (float(env['ALPHA'])+float(env["BUNDLEINDEX"])*float(env["DALPHA"])))) if "DALPHA" in env.keys() else env['ALPHA']
     npoly = str(int(env['NPOLY'])+int(env["BUNDLEINDEX"])*int(env["DNPOLY"])) if "DNPOLY" in env.keys() else env['NPOLY']
     nm = str(int(env['NM'])+int(env["BUNDLEINDEX"])*int(env["DNM"])) if "DNM" in env.keys() else env['NM']

     module = copy.deepcopy(moduleTemplates["IRL"])
     module["id"]["name"] = "epack"
     module["options"]["op"] = "stag_op_schur"
     module["options"]["lanczosParams"]["Cheby"]["alpha"] = '0'
     module["options"]["lanczosParams"]["Cheby"]["beta"] = '10' #env['BETA']
     module["options"]["lanczosParams"]["Cheby"]["Npoly"] = '11'
     module["options"]["lanczosParams"]["Nstop"] = '1000'#env['NSTOP']
     module["options"]["lanczosParams"]["Nk"] = '1050'#env['NK']
     module["options"]["lanczosParams"]["Nm"] = '1100'
     module["options"]["output"] = "eigen/eig{ens}nv2000er8_grid_{series}".format(ens=env['ENS'],series=env['SERIES'])
     module["options"]["multiFile"] = 'true'
     module["options"]["lanczosParams"]["resid"] = '1e-1'
                    
     modules.append(module)
    
     params["grid"]["modules"] = {"module":modules}

     return params
