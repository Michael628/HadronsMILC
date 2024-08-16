import copy
import os
import random

def buildParams(**moduleTemplates):

     env = os.environ

     time = int(env['TIME'])
     gammas = env["GAMMAS"].split()
     output = env["OUTPUT"].format(series=env['SERIES'],ens=env['ENS'],eigs=env['EIGS'],dt=env['DT'],noise=env['NOISE'])
     os.makedirs(output,exist_ok=True)
     params = {
          "grid":{
               "global":{
                    "nt":time,
                    "trajCounter":{
                         "start":env["CFG"],
                         "end":"10000",
                         "step":"10000",
                    },
                    "diskVectorDir":env["DISKVECTOR"],
                    "output":output,
               },
               "modules":{},
          },
     }

     a2aMatrix = []
     product = []
     Nh = int(env['NOISE'])*3*time
     
     doubleeigs=2*int(env['EIGS'])
     elabel=f"e{env['EIGS']}"

     seed=env["SEEDSTRING"]
     for i,sw in enumerate(range(int(env["NSEEDS"]))):
          wseed=f"w{seed}{sw}"

          for j,sv in enumerate(range(int(env["NSEEDS"]))):
               vseed=f"v{seed}{sv}"

               for gamma in gammas:

                    if j == 0 and i == 0:
                         elem = copy.deepcopy(moduleTemplates["contractA2Amatrix"])
                         elem["file"] = f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/proj_mesons/mf_{env['SERIES']}_{elabel}_{elabel}.{env['CFG']}/{gamma}_0_0_0.h5"
                         elem["dataset"] = f"{gamma}_0_0_0"
                         elem["cacheSize"] = "1"
                         elem["ni"] = doubleeigs
                         elem["nj"] = doubleeigs
                         elem["niOffset"] = 0
                         elem["njOffset"] = 0 
                         elem["name"] = gamma+elabel+elabel
                         a2aMatrix.append(elem)

                    if j == 0:
                         elem = copy.deepcopy(moduleTemplates["contractA2Amatrix"])
                         elem["file"] = f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/proj_mesons/mf_{env['SERIES']}_{wseed}_{elabel}.{env['CFG']}/{gamma}_0_0_0.h5"
                         elem["dataset"] = f"{gamma}_0_0_0"
                         elem["cacheSize"] = "1"
                         elem["ni"] = Nh
                         elem["nj"] = doubleeigs
                         elem["niOffset"] = 0
                         elem["njOffset"] = 0 
                         elem["name"] = gamma+wseed+elabel
                         a2aMatrix.append(elem)

                    if i == 0:
                         elem = copy.deepcopy(moduleTemplates["contractA2Amatrix"])
                         elem["file"] = f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/proj_mesons/mf_{env['SERIES']}_{elabel}_{vseed}.{env['CFG']}/{gamma}_0_0_0.h5"
                         elem["dataset"] = f"{gamma}_0_0_0"
                         elem["cacheSize"] = "1"
                         elem["ni"] = doubleeigs
                         elem["nj"] = Nh
                         elem["niOffset"] = 0
                         elem["njOffset"] = 0 
                         elem["name"] = gamma+elabel+vseed
                         a2aMatrix.append(elem)

                    if sw == sv:
                         continue

                    elem = copy.deepcopy(moduleTemplates["contractA2Amatrix"])
                    elem["file"] = f"e{env['EIGS']}n{env['NOISE']}dt{env['DT']}/proj_mesons/mf_{env['SERIES']}_{wseed}_{vseed}.{env['CFG']}/{gamma}_0_0_0.h5"
                    elem["dataset"] = f"{gamma}_0_0_0"
                    elem["cacheSize"] = "1"
                    elem["ni"] = Nh
                    elem["nj"] = Nh
                    elem["niOffset"] = 0
                    elem["njOffset"] = 0 
                    elem["name"] = gamma+wseed+vseed
                    a2aMatrix.append(elem)


     for i,s1 in enumerate(range(int(env["NSEEDS"]))):
          wseed1=f"w{seed}{s1}"
          vseed1=f"v{seed}{s1}"

          for j,s2 in enumerate(range(int(env["NSEEDS"]))):
               wseed2=f"w{seed}{s2}"
               vseed2=f"v{seed}{s2}"

               for gamma in gammas:
                    if j == 0:
                         elem = copy.deepcopy(moduleTemplates["contractA2Aproduct"])
                         elem["terms"] = f"{gamma+wseed1+elabel} {gamma+elabel+vseed1}"
                         elem["translations"] = f"0..{time-1}"
                         product.append(elem)

                    if i == 0:
                         elem = copy.deepcopy(moduleTemplates["contractA2Aproduct"])
                         elem["terms"] = f"{gamma+elabel+vseed2} {gamma+wseed2+elabel}"
                         elem["translations"] = f"0..{time-1}"
                         product.append(elem)

                    if s2 <= s1:
                         continue
                    
                    elem = copy.deepcopy(moduleTemplates["contractA2Aproduct"])
                    elem["terms"] = f"{gamma+wseed1+vseed2} {gamma+wseed2+vseed1}"
                    elem["translations"] = f"0..{time-1}"
                    product.append(elem)

     params["grid"]["a2aMatrix"] = {"elem":a2aMatrix}
     params["grid"]["product"] = {"elem":product}
     
     return params
